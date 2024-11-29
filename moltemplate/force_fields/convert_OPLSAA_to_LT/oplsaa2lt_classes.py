#!/usr/bin/env python3

from abc import ABC, abstractmethod
from oplsaa2lt_utils import data_from_atm_num


class Atom:
    """
    An Atom object, representing a single OPLSAA atom type,
    based on the paper Jorgensen_et_al-2024-The_Journal_of_Physical_Chemistry_B.
    NB: some values like charge-sigma-epsilon are kept as strings, to assure
    exactly the same values from the origianl FF file are taken

    Args and Attributes:
        type_id (int): the integer type ID from OPLSAA
        atomic_number (int): the atomic number of the atom
        type_str (str): the string representation of the atom type "subclass"
        charge (str): the charge of the atom
        sigma (str): the LJ-sigma value for the atom
        epsilon (str): the LJ-epsilon value for the atom
        comment (str): a comment associated with the atom (usually read from the original FF file)

    Attributes:
        element (str): the element of the atom (inferred from atomic_number)
        mass (str): the mass of the atom (inferred from atomic_number)

    Returns:
        an instance of the Atom class with the provided parameters
    """
    def __init__(self,
                 type_id: int,
                 atomic_number: int,
                 type_str: str,
                 charge: str,
                 sigma: str,
                 epsilon: str,
                 comment: str):

        self.type_id = type_id
        self.atomic_number = atomic_number
        self.type_str = type_str
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon
        self.bonded_type = type_str

        # NOTE: taking care of LP, which has atomic number 02, but it's not He...
        if self.type_str == "LP":
            self.atomic_number = 0

        self.element: str = data_from_atm_num[atomic_number]['element']
        self.mass: str = data_from_atm_num[atomic_number]['mass']

        self.comment = f"{self.element:2s} - {self.bonded_type:4s} | {comment}"

    @property
    def type_long(self) -> str:
        """
        This property returns a string representing the full type identifier for the atom.
        It includes the atom's type ID, followed by the bonded type for each type of bonded
        interaction (bond, angle, dihedral, improper), underscore separated.
        """
        l = f"{self.type_id}"
        l += f"_b{self.bonded_type}"
        l += f"_a{self.bonded_type}"
        l += f"_d{self.bonded_type}"
        l += f"_i{self.bonded_type}"
        return l

    @property
    def charge_line(self) -> str:
        return f"set type @atom:{self.type_id:<4d} charge {self.charge:>10s} # {self.comment}\n"

    @property
    def mass_line(self) -> str:
        return f"@atom:{self.type_id:<4d} {self.mass:<s}\n"

    @property
    def repl_line(self) -> str:
        return f"replace{{ @atom:{self.type_id:<4d} @atom:{self.type_long} }}\n"

    @property
    def nb_line(self) -> str:
        l = f"pair_coeff @atom:{self.type_long} @atom:{self.type_long}"
        l += f" {self.epsilon} {self.sigma}\n"
        return l


class BondedInteraction(ABC):
    """A base class for bonded interactions, providing common functionalities"""

    kind: str

    def __init__(self, types: list[str], comment: str):
        self.types = [ty.strip() for ty in types]
        self.comment = comment
        self.duplicate_count = 0
        self.to_skip = False

    @property
    def ty1(self) -> str:
        return self.types[0]

    @property
    def ty2(self) -> str:
        return self.types[1]

    @property
    def ty3(self) -> str:
        return self.types[2]

    @property
    def ty4(self) -> str:
        return self.types[3]

    @property
    def num_types(self) -> int:
        return len(self.types)

    @property
    def typename(self) -> str:
        """ This property returns a string with the name of the bonded interaction,
        composed by the names of the "base" atom types that form such interaction,
        separated by an underscore.

        NOTE: the bonded interaction name can't have "*" or "?" characters,
        so here such characters will be replaced with another character that,
        at the time of writing, is not used for atom types.
        """
        renamed_types = []
        for ty in self.types:
            renamed_types.append(ty.replace("*", "£").replace("?", "€"))
        return "_".join(renamed_types)

    @property
    def _coeff_line_base(self):
        return f"{type(self).kind}_coeff @{type(self).kind}:{self.typename}"

    @property
    @abstractmethod
    def coeff_line(self) -> str:
        pass

    @property
    @abstractmethod
    def bytype_line(self) -> str:
        pass

    @property
    def sort_key(self):
        """ This property returns a tuple that can be used to sort the interactions
        in a way that the more general interactions (the ones with * and ?)
        are placed at the beginning of the section."""

        def prioritize(ty):
            priority = lowest_priority
            if ty == "*":
                priority = 0
            elif ty == "?":
                priority = 1
            elif ty == "**":
                priority = 2
            elif ty == "??":
                priority = 3
            elif ty.startswith("*"):
                priority = 4
            elif ty.startswith("?"):
                priority = 5
            elif ty.endswith("*"):
                priority = 6
            elif ty.endswith("?"):
                priority = 7
            return priority

        lowest_priority = 8
        priority_tuple = (
            prioritize(self.ty1),
            prioritize(self.ty2),
            prioritize(self.ty3) if self.num_types >= 3 else lowest_priority,
            prioritize(self.ty4) if self.num_types >= 4 else lowest_priority
            )
        return priority_tuple

    @property
    def _check_if_comment(self) -> str:
        # If there are multiple verions of this interaction, comment out the duplicates
        if self.duplicate_count > 1:
            return "#"
        # The next 2 lines are no longer needed
        # (Instead, we replace "C(O)" with "CparenO")
        # if self.typename == "HC_CT_CT_C(O)"
        #     return "#"
        return ""

    @property
    def _duplicate_count_str(self) -> str:
        # If there N multiple verions of this interaction
        # acting on the same atom types, then
        # add a unique integer to the name of this interaction type.
        if self.duplicate_count > 0:
            return f"__{self.duplicate_count}"
        return ""


class Bond(BondedInteraction):
    kind = "bond"
    def __init__(self, types: list[str], k: str, eq: str, comment: str):
        super().__init__(types, comment)
        if len(self.types) != 2:
            raise ValueError(f"interaction of type '{type(self).kind}' must have 2 types...")
        self.k = k
        self.eq = eq
        self.comment = comment

    @property
    def bytype_line(self) -> str:
        l = f"{self._check_if_comment}"
        l += f"@{type(self).kind}:{self.typename}"
        l += f"{self._duplicate_count_str}"
        l += f" @atom:*_b{self.ty1}_a*_d*_i*"
        l += f" @atom:*_b{self.ty2}_a*_d*_i*\n"
        return l

    @property
    def coeff_line(self) -> str:
        l = ""  # l = f"{self._check_if_comment}"
        l += f"{self._coeff_line_base}{self._duplicate_count_str} {self.k} {self.eq} # {self.comment}\n"
        return l


class Angle(BondedInteraction):
    kind = "angle"
    def __init__(self, types: list[str], k: str, eq: str, comment: str):
        super().__init__(types, comment)
        if len(self.types) != 3:
            raise ValueError(f"interaction of type '{type(self).kind}' must have 3 types...")
        self.k = k
        self.eq = eq
        self.comment = comment

    @property
    def bytype_line(self) -> str:
        l = f"{self._check_if_comment}"
        l += f"@{type(self).kind}:{self.typename}"
        l += f"{self._duplicate_count_str}"
        l += f" @atom:*_b*_a{self.ty1}_d*_i*"
        l += f" @atom:*_b*_a{self.ty2}_d*_i*"
        l += f" @atom:*_b*_a{self.ty3}_d*_i*\n"
        return l

    @property
    def coeff_line(self) -> str:
        l = ""  # l = f"{self._check_if_comment}"
        l += f"{self._coeff_line_base}{self._duplicate_count_str} {self.k} {self.eq} # {self.comment}\n"
        return l


class Dihedral(BondedInteraction):
    kind = "dihedral"
    def __init__(self, types: list[str], v1, v2, v3, v4, comment: str):
        super().__init__(types, comment)
        if len(self.types) != 4:
            raise ValueError(f"interaction of type '{type(self).kind}' must have 4 types...")
        self.v1, self.v2, self.v3, self.v4 = v1, v2, v3, v4
        self.comment = comment

    @property
    def bytype_line(self) -> str:
        l = f"{self._check_if_comment}"
        l += f"@{type(self).kind}:{self.typename}"
        l += f"{self._duplicate_count_str}"
        l += f" @atom:*_b*_a*_d{self.ty1}_i*"
        l += f" @atom:*_b*_a*_d{self.ty2}_i*"
        l += f" @atom:*_b*_a*_d{self.ty3}_i*"
        l += f" @atom:*_b*_a*_d{self.ty4}_i*\n"
        return l

    @property
    def coeff_line(self) -> str:
        l = ""  # l = f"{self._check_if_comment}"
        l += f"{self._coeff_line_base}{self._duplicate_count_str} {self.v1} {self.v2} {self.v3} {self.v4}"
        l += f" # {self.comment} \n"
        return l


class Improper(BondedInteraction):
    kind = "improper"
    def __init__(self, types: list[str], v1, v2, v3, v4, comment: str):
        super().__init__(types, comment)
        if len(self.types) != 4:
            raise ValueError(f"interaction of type '{type(self).kind}' must have 4 types...")
        self.v1, self.v2, self.v3, self.v4 = v1, v2, v3, v4
        self.comment = comment

        # replacing every X/Y/Z with "*", I don't know if that's the correct approach,
        # as if this was the intended behaviour I think in the FF file they would have used
        # "*", as they did in the dihedrals sections...
        for i, ty in enumerate(self.types):
            if ty in ("X", "Y", "Z"):
                self.types[i] = "*"

    @property
    def bytype_line(self) -> str:
        l = f"{self._check_if_comment}"
        l += f"@{type(self).kind}:{self.typename}"
        l += f" @atom:*_b*_a*_d*_i{self.ty1}"
        l += f" @atom:*_b*_a*_d*_i{self.ty2}"
        l += f" @atom:*_b*_a*_d*_i{self.ty3}"
        l += f" @atom:*_b*_a*_d*_i{self.ty4}\n"
        return l

    @property
    def coeff_line(self) -> str:
        l = ""  # l = f"{self._check_if_comment}"
        # If using "improper_style cvff", then use:
        l += f"{self._coeff_line_base}{self._duplicate_count_str} {float(self.v2)/2:.4f} -1 2 # {self.comment}\n"
        # If using "improper_style harmonic", then use this instead:
        # l += f"{self._coeff_line_base}{self._duplicate_count_str} {float(self.v2)} 180.0 # {self.comment}\n"
        return l
