"""
mcr_ptm_annotator.ptm_database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Experimentally characterised post-translational modifications (PTMs)
in the McrA subunit of methyl-coenzyme M reductase.

All entries are sourced from published experimental data only.
No predicted or computationally inferred modifications are included.

Primary structural reference
-----------------------------
Ermler et al. (1997) Science 278:1457-1462
  Crystal structure of MCR from Methanothermobacter marburgensis, PDB: 1MRO.
  This is the reference numbering used for all positions below.
  Positions in other organisms will differ due to sequence variation.

PTM characterisation references
--------------------------------
Ermler et al. (1997) Science 278:1457      -- initial identification of 5 PTMs in PDB 1MRO
Selmer et al. (2000) J Biol Chem 275:3755  -- methylation cross-organism survey
Kahnt et al. (2007) FEBS J 274:4913        -- mass spec confirmation, cross-species survey
Nayak et al. (2017) eLife 6:e29218         -- genetic proof of thioglycine (tfuA/ycaO)
Wagner et al. (2016) Angew Chem 55:10630   -- didehydroaspartate, 6th PTM
Deobald et al. (2018) Sci Rep 8:7404       -- radical SAM methylase (mmpX) for Arg271
Nayak et al. (2020) PLoS Biol 18:e3000507  -- thioamidation required for MCR activity

Known variability across organisms (Kahnt 2007, Selmer 2000):
  - Thioglycine and 1-N-methylhistidine: conserved in all methanogens examined
  - 5-(S)-methylarginine: present in all methanogens examined, absent in ANME-1
  - 2-(S)-methylglutamine: absent in Methanosarcina barkeri (Selmer 2000)
  - S-methylcysteine: low abundance or absent in many methanogens, e.g. M. maripaludis
  - Didehydroaspartate: present in M. marburgensis and M. barkeri, absent in M. wolfeii
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class KnownPTM:
    """
    A single experimentally characterised PTM in McrA.

    Attributes
    ----------
    name : str
        Common name of the modification.
    residue_type : str
        Single-letter amino acid code of the unmodified residue.
    position_marburgensis : int
        1-based residue position in M. marburgensis McrA (PDB 1MRO).
        Use this ONLY as a reference point; positions in other organisms
        differ and must be mapped by sequence alignment.
    modification : str
        Chemical description of the modification.
    pdb_code : str
        Three-letter PDB HETATM code for the modified residue in 1MRO,
        or empty string if not in 1MRO.
    references : tuple[str, ...]
        Key citations for experimental characterisation.
    notes : str
        Additional context including known variability across organisms.
    """
    name: str
    residue_type: str
    position_marburgensis: int
    modification: str
    pdb_code: str
    references: tuple
    notes: str = ""


# ── Verified MCR PTMs (McrA subunit, M. marburgensis numbering) ────────────
#
# The first 5 modifications were identified in PDB 1MRO (Ermler 1997) and
# confirmed by mass spectrometry (Kahnt 2007).
# A 6th modification (didehydroaspartate) was identified in M. marburgensis
# by Wagner et al. (2016) using mass spectrometry and high-resolution X-ray
# crystallography (PDB 5A0Y). It is adjacent to thioglycine.
#
KNOWN_PTMS: tuple[KnownPTM, ...] = (

    KnownPTM(
        name="thioglycine",
        residue_type="G",
        position_marburgensis=445,
        modification="Thioamidation: backbone carbonyl oxygen replaced by sulfur",
        pdb_code="GL3",
        references=(
            "Ermler et al. (1997) Science 278:1457",
            "Nayak et al. (2017) eLife 6:e29218",
            "Nayak et al. (2020) PLoS Biol 18:e3000507",
        ),
        notes=(
            "Installed by YcaO/TfuA (ycaO-tfuA locus, Nayak 2017). "
            "In M. acetivorans this corresponds to Gly465. "
            "Loss causes severe growth defects on low-energy substrates "
            "and at elevated temperature (39-45°C, Nayak 2017). "
            "Conserved in all methanogens examined (Kahnt 2007)."
        ),
    ),

    KnownPTM(
        name="didehydroaspartate",
        residue_type="D",
        position_marburgensis=446,
        modification="Dehydration of aspartate: alpha-beta unsaturated amino acid",
        pdb_code="",  # Not in PDB 1MRO; identified in PDB 5A0Y (Wagner 2016)
        references=(
            "Wagner et al. (2016) Angew Chem Int Ed Engl 55:10630",
            "Nayak et al. (2017) eLife 6:e29218",
        ),
        notes=(
            "Adjacent to thioglycine (Gly445). Identified in M. marburgensis "
            "MCR I and II and in M. barkeri by mass spectrometry and X-ray "
            "crystallography (PDB 5A0Y, Wagner 2016). "
            "Absent in M. wolfeii — dispensable but may fine-tune catalytic "
            "efficiency (Wagner 2016). "
            "In M. acetivorans, corresponds to Asp470 (Nayak 2017). "
            "NOT present in PDB 1MRO; use PDB 5A0Y as reference for this PTM."
        ),
    ),

    KnownPTM(
        name="1-N-methylhistidine",
        residue_type="H",
        position_marburgensis=257,
        modification="Methylation of the N1 nitrogen of the histidine imidazole ring",
        pdb_code="MHS",
        references=(
            "Ermler et al. (1997) Science 278:1457",
            "Kahnt et al. (2007) FEBS J 274:4913",
        ),
        notes=(
            "Conserved in all methanogens examined (Kahnt 2007). "
            "Biosynthetic enzyme not fully characterised as of 2024. "
            "Proposed to alter pKa and position the imidazole ring "
            "for coenzyme B binding (Grabarse et al. 2000)."
        ),
    ),

    KnownPTM(
        name="5-(S)-methylarginine",
        residue_type="R",
        position_marburgensis=271,
        modification="sp3-C methylation at the C-delta of arginine",
        pdb_code="AGM",
        references=(
            "Ermler et al. (1997) Science 278:1457",
            "Kahnt et al. (2007) FEBS J 274:4913",
            "Deobald et al. (2018) Sci Rep 8:7404",
        ),
        notes=(
            "Installed by the radical SAM methyltransferase encoded by mmpX "
            "(Deobald 2018). Present in all methanogens examined, absent in "
            "ANME-1 (Kahnt 2007). In M. maripaludis, the equivalent residue "
            "is Arg275; loss of methylation profoundly reduces methanogenesis "
            "and growth (Lyu et al. 2020 J Bacteriol 202:e00654-19)."
        ),
    ),

    KnownPTM(
        name="2-(S)-methylglutamine",
        residue_type="Q",
        position_marburgensis=400,
        modification="sp3-C methylation at the C-alpha of glutamine",
        pdb_code="MGN",
        references=(
            "Ermler et al. (1997) Science 278:1457",
            "Kahnt et al. (2007) FEBS J 274:4913",
            "Selmer et al. (2000) J Biol Chem 275:3755",
        ),
        notes=(
            "Absent in Methanosarcina barkeri (Selmer 2000, Lyu et al. 2018). "
            "Biosynthetic enzyme not fully characterised as of 2024."
        ),
    ),

    KnownPTM(
        name="S-methylcysteine",
        residue_type="C",
        position_marburgensis=452,
        modification="S-methylation of the cysteine thiol",
        pdb_code="SMC",
        references=(
            "Ermler et al. (1997) Science 278:1457",
            "Kahnt et al. (2007) FEBS J 274:4913",
            "Selmer et al. (2000) J Biol Chem 275:3755",
        ),
        notes=(
            "Variable: low abundance or absent in many methanogens, "
            "including M. maripaludis (Kahnt 2007). "
            "Present in M. marburgensis (PDB 1MRO). "
            "ANME-1 MCR shows a different PTM pattern at this position "
            "(Shima et al. 2012 Nature 481:98)."
        ),
    ),
)
