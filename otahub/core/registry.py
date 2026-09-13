"""Discovery and lookup across every archetype spec on disk."""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Iterator

from .archetype import Archetype
from .spec import ArchetypeSpec, FamilySpec


def default_spec_dir() -> Path:
    """Where specs live. Overridable with ``OTAHUB_SPEC_DIR`` for tests."""
    env = os.environ.get("OTAHUB_SPEC_DIR")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / "specs"


class Registry:
    """An in-memory catalogue of antenna archetypes loaded from JSON specs."""

    def __init__(self, archetypes: Iterable[Archetype] = (),
                 load_errors: dict[str, str] | None = None) -> None:
        self._by_key: dict[str, Archetype] = {}
        self.load_errors: dict[str, str] = dict(load_errors or {})
        for a in archetypes:
            self.add(a)

    # ------------------------------------------------------------- building

    def add(self, archetype: Archetype) -> None:
        if archetype.key in self._by_key:
            existing = self._by_key[archetype.key].spec.source_file
            raise ValueError(
                f"duplicate archetype key {archetype.key!r} "
                f"(already defined in {existing})"
            )
        self._by_key[archetype.key] = archetype

    @classmethod
    def load(cls, spec_dir: str | Path | None = None) -> "Registry":
        """Load every ``*.json`` in `spec_dir`.

        A file that fails to parse is recorded in ``load_errors`` rather than
        taking the whole catalogue down — one bad spec should not make the
        tool unusable.
        """
        directory = Path(spec_dir) if spec_dir else default_spec_dir()
        registry = cls()
        if not directory.is_dir():
            return registry
        for path in sorted(directory.glob("*.json")):
            try:
                family = FamilySpec.load(path)
            except Exception as exc:  # noqa: BLE001 - surfaced to the user
                registry.load_errors[path.name] = f"{type(exc).__name__}: {exc}"
                continue
            for spec in family:
                try:
                    registry.add(Archetype(spec))
                except ValueError as exc:
                    registry.load_errors[f"{path.name}:{spec.key}"] = str(exc)
        return registry

    # -------------------------------------------------------------- querying

    def __len__(self) -> int: return len(self._by_key)
    def __iter__(self) -> Iterator[Archetype]: return iter(self._by_key.values())
    def __contains__(self, key: object) -> bool: return key in self._by_key

    def __getitem__(self, key: str) -> Archetype:
        try:
            return self._by_key[key]
        except KeyError:
            near = self.suggest(key)
            hint = f" Did you mean: {', '.join(near)}?" if near else ""
            raise KeyError(f"no archetype named {key!r}.{hint}") from None

    def get(self, key: str, default: Any = None) -> Archetype | Any:
        return self._by_key.get(key, default)

    @property
    def keys(self) -> list[str]:
        return sorted(self._by_key)

    @property
    def families(self) -> list[str]:
        return sorted({a.family for a in self})

    def by_family(self, family: str) -> list[Archetype]:
        return sorted((a for a in self if a.family == family), key=lambda a: a.key)

    def covering(self, f_hz: float) -> list[Archetype]:
        """Archetypes whose stated validity band contains `f_hz`."""
        return sorted((a for a in self if a.covers(f_hz)), key=lambda a: a.key)

    def search(self, text: str) -> list[Archetype]:
        """Free-text match over key, name and summary."""
        needle = text.lower().strip()
        if not needle:
            return sorted(self, key=lambda a: a.key)
        hits = []
        for a in self:
            haystack = f"{a.key} {a.name} {a.spec.summary} {a.family}".lower()
            if needle in haystack:
                hits.append(a)
        return sorted(hits, key=lambda a: a.key)

    def suggest(self, key: str, n: int = 3) -> list[str]:
        import difflib
        return difflib.get_close_matches(key, self.keys, n=n, cutoff=0.5)

    # ------------------------------------------------------------ validation

    def problems(self) -> dict[str, list[str]]:
        """Static faults per archetype, plus any file-level load errors."""
        out: dict[str, list[str]] = {}
        for name, err in self.load_errors.items():
            out.setdefault(f"<file {name}>", []).append(err)
        for a in self:
            probs = a.spec.problems()
            if probs:
                out[a.key] = probs
        return out


@lru_cache(maxsize=1)
def default_registry() -> Registry:
    """Process-wide registry loaded from the default spec directory."""
    return Registry.load()
