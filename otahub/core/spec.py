"""Typed model of the archetype spec files under ``specs/``.

Loading is deliberately *tolerant about shape* and *strict about meaning*:
optional fields may be absent, but a malformed formula or an unresolvable
dependency is collected as a hard problem rather than silently dropped.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Sequence

from .expr import ExprError, referenced_symbols


@dataclass(frozen=True, slots=True)
class ParameterSpec:
    symbol: str
    name: str = ""
    unit: str = ""
    role: str = "geometry"   # geometry | requirement | material | assumption | derived
    description: str = ""
    typical: str = ""

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "ParameterSpec":
        return cls(
            symbol=d["symbol"],
            name=d.get("name", d["symbol"]),
            unit=d.get("unit", ""),
            role=d.get("role", "geometry"),
            description=d.get("description", ""),
            typical=str(d.get("typical", "")),
        )


@dataclass(frozen=True, slots=True)
class SynthesisRule:
    output: str
    expr: str
    units_out: str = ""
    notes: str = ""
    depends_on: tuple[str, ...] = ()

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "SynthesisRule":
        return cls(
            output=d["output"], expr=d["expr"],
            units_out=d.get("units_out", ""), notes=d.get("notes", ""),
            depends_on=tuple(d.get("depends_on", ())),
        )


@dataclass(frozen=True, slots=True)
class AnalysisRule:
    metric: str
    expr: str
    units_out: str = ""
    notes: str = ""
    depends_on: tuple[str, ...] = ()

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "AnalysisRule":
        return cls(
            metric=d["metric"], expr=d["expr"],
            units_out=d.get("units_out", ""), notes=d.get("notes", ""),
            depends_on=tuple(d.get("depends_on", ())),
        )


@dataclass(frozen=True, slots=True)
class KnownCase:
    """A citable numeric fact that becomes a pytest case."""
    given: dict[str, float]
    expect: dict[str, float]
    tol_pct: float = 5.0
    source: str = ""
    #: Absolute tolerance, used INSTEAD of tol_pct when supplied. Required for
    #: any expectation of zero, where relative error is undefined - a resonant
    #: reactance of 3 ohm against an expected 0 is excellent, not a 300% miss.
    tol_abs: float | None = None

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "KnownCase":
        tol_abs = d.get("tol_abs")
        return cls(
            given=dict(d.get("given", {})),
            expect=dict(d.get("expect", {})),
            tol_pct=float(d.get("tol_pct", 5.0)),
            source=d.get("source", ""),
            tol_abs=None if tol_abs is None else float(tol_abs),
        )


@dataclass(slots=True)
class ArchetypeSpec:
    key: str
    name: str
    family: str
    summary: str = ""
    freq_range_hz: tuple[float, float] = (0.0, float("inf"))
    parameters: tuple[ParameterSpec, ...] = ()
    synthesis: tuple[SynthesisRule, ...] = ()
    analysis: tuple[AnalysisRule, ...] = ()
    validity: tuple[str, ...] = ()
    known_cases: tuple[KnownCase, ...] = ()
    references: tuple[str, ...] = ()
    confidence: str = "normal"
    tables: dict[str, Any] = field(default_factory=dict)
    source_file: str = ""

    @property
    def symbols(self) -> frozenset[str]:
        """Every symbol this archetype may legally reference."""
        return frozenset(
            [p.symbol for p in self.parameters]
            + [r.output for r in self.synthesis]
            + [r.metric for r in self.analysis]
            + ["f0", "lambda0", "k0"]
        )

    @classmethod
    def from_json(cls, d: dict[str, Any], family: str, source_file: str = "") -> "ArchetypeSpec":
        rng = d.get("freq_range_hz") or [0.0, float("inf")]
        return cls(
            key=d["key"],
            name=d.get("name", d["key"].replace("_", " ").title()),
            family=d.get("family", family),
            summary=d.get("summary", ""),
            freq_range_hz=(float(rng[0]), float(rng[1])),
            parameters=tuple(ParameterSpec.from_json(p) for p in d.get("parameters", [])),
            synthesis=tuple(SynthesisRule.from_json(r) for r in d.get("synthesis", [])),
            analysis=tuple(AnalysisRule.from_json(r) for r in d.get("analysis", [])),
            validity=tuple(d.get("validity", [])),
            known_cases=tuple(KnownCase.from_json(k) for k in d.get("known_cases", [])),
            references=tuple(d.get("references", [])),
            confidence=d.get("confidence", "normal"),
            tables=dict(d.get("tables", {})),
            source_file=source_file,
        )

    # ------------------------------------------------------------- validation

    def problems(self) -> list[str]:
        """Static faults in this archetype. Empty list means structurally sound."""
        out: list[str] = []
        syms = self.symbols
        seen_outputs: set[str] = set()

        for rule in self.synthesis:
            if rule.output in seen_outputs:
                out.append(f"synthesis output {rule.output!r} defined more than once")
            seen_outputs.add(rule.output)
            try:
                used = referenced_symbols(rule.expr, syms)
            except ExprError as exc:
                out.append(f"synthesis {rule.output!r}: {exc}")
                continue
            unknown = used - syms
            if unknown:
                out.append(f"synthesis {rule.output!r} references unknown {sorted(unknown)}")
            if rule.depends_on:
                declared = set(rule.depends_on)
                missing = used - declared
                if missing:
                    out.append(
                        f"synthesis {rule.output!r} reads {sorted(missing)} "
                        f"but does not declare them in depends_on"
                    )

        for rule in self.analysis:
            try:
                used = referenced_symbols(rule.expr, syms)
            except ExprError as exc:
                out.append(f"analysis {rule.metric!r}: {exc}")
                continue
            unknown = used - syms
            if unknown:
                out.append(f"analysis {rule.metric!r} references unknown {sorted(unknown)}")

        if not self.known_cases:
            out.append("no known_cases: nothing about this archetype is verifiable")
        for i, case in enumerate(self.known_cases):
            if not case.expect:
                out.append(f"known case {i} expects nothing: it would pass untested")
            for name, want in case.expect.items():
                if isinstance(want, bool) or not isinstance(want, (int, float, complex)):
                    out.append(f"known case {i}: {name!r} expects {want!r}, not a number")
                elif want == 0 and case.tol_abs is None:
                    out.append(f"known case {i}: {name!r} expects zero without tol_abs; "
                               f"relative error has no meaning there")
        if self.freq_range_hz[0] >= self.freq_range_hz[1]:
            out.append(f"degenerate freq_range_hz {self.freq_range_hz}")
        return out


@dataclass(slots=True)
class FamilySpec:
    family: str
    archetypes: tuple[ArchetypeSpec, ...]
    source_file: str = ""

    def __iter__(self) -> Iterator[ArchetypeSpec]:
        return iter(self.archetypes)

    def __len__(self) -> int:
        return len(self.archetypes)

    @classmethod
    def load(cls, path: str | Path) -> "FamilySpec":
        path = Path(path)
        raw = json.loads(path.read_text())
        family = raw.get("family", path.stem)
        items = raw.get("archetypes", [])
        return cls(
            family=family,
            archetypes=tuple(
                ArchetypeSpec.from_json(a, family, str(path)) for a in items
            ),
            source_file=str(path),
        )
