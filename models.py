from dataclasses import dataclass

@dataclass(frozen=True)
class Species:
    smiles: str
    experimental_hof: float
    uncertainty: float
    calculated_hof: float

@dataclass(frozen=True)
class Reaction:
    reactants: tuple[tuple[int, str], ...]
    products: tuple[tuple[int, str], ...]

@dataclass(frozen=True)
class FeatureOptions:
    electron_pair: bool = True
    lone_electron: bool = True
    hydro_bond: bool = True
    normal_bond: bool = False