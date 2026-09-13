"""The runtime object: a spec plus the ability to synthesise and analyse."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from .constants import C0
from .expr import ExprError, evaluate, referenced_symbols
from .spec import ArchetypeSpec


class SynthesisError(RuntimeError):
    """Raised when requirements are insufficient or inconsistent."""


@dataclass(slots=True)
class DesignResult:
    """Outcome of a synthesis run, with enough provenance to audit it."""
    archetype: str
    requirements: dict[str, Any]
    parameters: dict[str, Any]
    metrics: dict[str, Any]
    units: dict[str, str]
    warnings: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    #: output symbol -> the requirement symbols that would unlock it
    unresolved: dict[str, list[str]] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.unresolved is None:
            self.unresolved = {}

    def missing_requirements(self) -> list[str]:
        """Every symbol that, if supplied, would unlock something unresolved."""
        seen: set[str] = set()
        for needs in (self.unresolved or {}).values():
            seen.update(needs)
        return sorted(seen - set(self.parameters) - set(self.requirements))

    def get(self, name: str, default: Any = None) -> Any:
        if name in self.parameters:
            return self.parameters[name]
        return self.metrics.get(name, default)

    def as_dict(self) -> dict[str, Any]:
        return {
            "archetype": self.archetype,
            "requirements": self.requirements,
            "parameters": self.parameters,
            "metrics": self.metrics,
            "units": self.units,
            "warnings": list(self.warnings),
        }


class Archetype:
    """Wraps an :class:`ArchetypeSpec` with the synthesis/analysis engine."""

    def __init__(self, spec: ArchetypeSpec) -> None:
        self.spec = spec

    # ------------------------------------------------------------- properties
    @property
    def key(self) -> str: return self.spec.key
    @property
    def name(self) -> str: return self.spec.name
    @property
    def family(self) -> str: return self.spec.family

    def __repr__(self) -> str:
        return f"<Archetype {self.spec.key} ({self.spec.family})>"

    def covers(self, f_hz: float) -> bool:
        lo, hi = self.spec.freq_range_hz
        return lo <= f_hz <= hi

    # ------------------------------------------------------------- synthesis

    def _seed(self, requirements: Mapping[str, Any]) -> dict[str, Any]:
        """Requirements plus the universally-derived frequency quantities."""
        known: dict[str, Any] = dict(requirements)
        # Material properties have defensible defaults (copper, mu_r = 1), as
        # do parameters explicitly declared an "assumption" (a thin-wire
        # radius, say). A design *requirement* has no defensible default, so
        # it is never silently invented - that would change the antenna.
        for p in self.spec.parameters:
            if p.role not in ("material", "assumption") or p.symbol in known:
                continue
            try:
                known[p.symbol] = float(p.typical)
            except (TypeError, ValueError):
                continue
        f0 = known.get("f0")
        if f0 is not None and f0 > 0:
            known.setdefault("lambda0", C0 / f0)
            known.setdefault("k0", 2.0 * math.pi * f0 / C0)
        return known

    def synthesize(self, strict: bool = False, **requirements: Any) -> DesignResult:
        """Requirements -> physical dimensions, then analyse the result.

        Rules are resolved by repeated relaxation on their *actual* referenced
        symbols rather than the declared ``depends_on``, so a spec with a stale
        dependency list still resolves correctly (the discrepancy is reported
        by :meth:`ArchetypeSpec.problems`, not silently honoured here).

        Resolution is *partial by design*. An archetype commonly carries
        optional rules — a wire radius derived from an efficiency target, a
        tuning capacitor — that need requirements a given caller has no reason
        to supply. Those become warnings and populate
        :attr:`DesignResult.unresolved`; everything that *can* be computed
        still is. Pass ``strict=True`` to demand a complete geometry instead.
        """
        known = self._seed(requirements)
        warnings: list[str] = []
        syms = self.spec.symbols | frozenset(known)

        pending = list(self.spec.synthesis)
        units: dict[str, str] = {p.symbol: p.unit for p in self.spec.parameters}

        progressed = True
        while pending and progressed:
            progressed = False
            still: list[Any] = []
            for rule in pending:
                try:
                    needed = referenced_symbols(rule.expr, syms)
                except ExprError as exc:
                    warnings.append(f"skipped {rule.output!r}: {exc}")
                    continue
                if needed <= set(known):
                    try:
                        known[rule.output] = evaluate(rule.expr, known)
                    except ExprError as exc:
                        warnings.append(f"{rule.output!r} failed: {exc}")
                        continue
                    if rule.units_out:
                        units[rule.output] = rule.units_out
                    progressed = True
                else:
                    still.append(rule)
            pending = still

        unresolved: dict[str, list[str]] = {}
        for rule in pending:
            missing = sorted(referenced_symbols(rule.expr, syms) - set(known))
            unresolved[rule.output] = missing
            warnings.append(
                f"{rule.output!r} not computed: needs {missing}. "
                f"Supply them as requirements to complete the design."
            )
        if strict and unresolved:
            detail = "; ".join(f"{k} needs {v}" for k, v in unresolved.items())
            raise SynthesisError(f"{self.spec.key}: incomplete design — {detail}")

        params = {k: v for k, v in known.items() if k not in requirements}
        metrics, mwarn = self._analyse_into(known, units, unresolved)
        warnings.extend(mwarn)
        warnings.extend(self._validity_warnings(known))

        return DesignResult(
            archetype=self.spec.key,
            requirements=dict(requirements),
            parameters=params,
            metrics=metrics,
            units=units,
            warnings=tuple(warnings),
            notes=tuple(self.spec.validity),
            unresolved=dict(unresolved),
        )

    # -------------------------------------------------------------- analysis

    def _analyse_into(
        self, known: Mapping[str, Any], units: dict[str, str],
        unresolved: dict[str, list[str]] | None = None,
    ) -> tuple[dict[str, Any], list[str]]:
        metrics: dict[str, Any] = {}
        warnings: list[str] = []
        scope = dict(known)
        syms = self.spec.symbols | frozenset(scope)

        # Metrics may depend on one another; relax until stable.
        pending = list(self.spec.analysis)
        progressed = True
        while pending and progressed:
            progressed = False
            still = []
            for rule in pending:
                try:
                    needed = referenced_symbols(rule.expr, syms)
                except ExprError as exc:
                    warnings.append(f"analysis {rule.metric!r}: {exc}")
                    continue
                if needed <= set(scope):
                    try:
                        value = evaluate(rule.expr, scope)
                    except ExprError as exc:
                        warnings.append(f"analysis {rule.metric!r}: {exc}")
                        continue
                    metrics[rule.metric] = value
                    scope[rule.metric] = value
                    if rule.units_out:
                        units[rule.metric] = rule.units_out
                    progressed = True
                else:
                    still.append(rule)
            pending = still

        for rule in pending:
            missing = sorted(referenced_symbols(rule.expr, syms) - set(scope))
            warnings.append(f"analysis {rule.metric!r} unavailable, needs {missing}")
            if unresolved is not None:
                unresolved[rule.metric] = missing
        return metrics, warnings

    def analyze(self, **params: Any) -> dict[str, Any]:
        """Analyse a hand-supplied geometry without running synthesis."""
        known = self._seed(params)
        units: dict[str, str] = {}
        metrics, _ = self._analyse_into(known, units)
        return metrics

    # -------------------------------------------------------------- checking

    def _validity_warnings(self, known: Mapping[str, Any]) -> list[str]:
        out: list[str] = []
        f0 = known.get("f0")
        if f0 is not None:
            lo, hi = self.spec.freq_range_hz
            if not (lo <= f0 <= hi):
                out.append(
                    f"f0={f0:.4g} Hz is outside this archetype's stated validity "
                    f"band [{lo:.4g}, {hi:.4g}] Hz"
                )
        if self.spec.confidence == "low":
            out.append(
                "spec is marked low-confidence: treat these numbers as indicative "
                "and verify in a full-wave solver before committing"
            )
        return out

    def check_known_cases(self) -> list[dict[str, Any]]:
        """Run every citable known case. Returns one record per expectation."""
        results: list[dict[str, Any]] = []
        for case in self.spec.known_cases:
            try:
                design = self.synthesize(**case.given)
            except Exception as exc:  # noqa: BLE001 - reported, not raised
                for name, want in case.expect.items():
                    results.append({
                        "archetype": self.spec.key, "quantity": name,
                        "expected": want, "actual": None, "error_pct": None,
                        "tol_pct": case.tol_pct, "tol_abs": case.tol_abs, "passed": False,
                        "source": case.source, "detail": f"{type(exc).__name__}: {exc}",
                    })
                continue
            for name, want in case.expect.items():
                got = _lookup(design, name)
                passed, err = _compare(got, want, case.tol_pct, case.tol_abs)
                results.append({
                    "archetype": self.spec.key, "quantity": name,
                    "expected": want, "actual": got, "error_pct": err,
                    "tol_pct": case.tol_pct, "tol_abs": case.tol_abs, "passed": passed,
                    "source": case.source,
                    "detail": "" if got is not None else "quantity not produced",
                })
        return results


# ------------------------------------------------------------------ helpers

_UNIT_SUFFIXES = ("_m", "_mm", "_cm", "_hz", "_ghz", "_mhz", "_ohm",
                  "_dbi", "_db", "_deg", "_rad", "_pct", "_percent")


def _lookup(design: DesignResult, name: str) -> Any:
    """Resolve an expectation key against a design result, tolerating suffixes.

    Specs write expectations as ``L_m`` or ``directivity_dbi``; the engine may
    hold them as ``L`` or ``directivity_dbi``. Try the literal name first, then
    progressively strip a trailing unit suffix.
    """
    value = design.get(name)
    if value is not None:
        return value

    lowered = name.lower()
    for suffix in _UNIT_SUFFIXES:
        if lowered.endswith(suffix):
            value = design.get(name[: -len(suffix)])
            if value is not None:
                return value

    # Fall back to stripping trailing underscore-separated segments, which
    # catches ad-hoc unit tags a spec may carry ("A_m2", "ka_equiv_-").
    # Only a candidate that actually resolves is accepted, so this widens
    # what matches without inventing a value.
    parts = name.split("_")
    for cut in (1, 2):
        if len(parts) > cut:
            value = design.get("_".join(parts[:-cut]))
            if value is not None:
                return value
    return None


def _compare(got: Any, want: Any, tol_pct: float,
             tol_abs: float | None = None) -> tuple[bool, float | None]:
    """Compare a computed value against an expectation.

    Uses absolute tolerance when the spec supplies one, which is the only
    meaningful test against an expected zero: relative error has no definition
    there, and a 3 ohm residual reactance where zero was wanted is a good
    result, not a 300% failure.
    """
    if got is None:
        return False, None
    try:
        if isinstance(got, complex) or isinstance(want, complex):
            g, w = complex(got), complex(want)
            delta = abs(g - w)
            denom = abs(w) if abs(w) > 0 else 1.0
        else:
            g, w = float(got), float(want)
            if not (math.isfinite(g) and math.isfinite(w)):
                return g == w, None
            delta = abs(g - w)
            denom = abs(w) if abs(w) > 1e-15 else 1.0
        err = delta / denom * 100.0
    except (TypeError, ValueError):
        return got == want, None
    if tol_abs is not None:
        return delta <= tol_abs, err
    return err <= tol_pct, err
