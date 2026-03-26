# Protein statistics — high-performance implementations of ExPASy ProtParam-style analyses
# All lookup tables are 256-element byte-indexed arrays for O(1) access (no Dict overhead).
# Optional GPU support is intentionally kept off the startup path.

# ─── Amino acid mass tables (indexed by UInt8 of uppercase ASCII letter) ───────

const _AA_MASS_MONO = fill(0.0, 256)
const _AA_MASS_AVG  = fill(0.0, 256)

# Monoisotopic masses (Da)
for (aa, mass) in [
    ('A', 71.03711), ('R', 156.10111), ('N', 114.04293), ('D', 115.02694),
    ('C', 103.00919), ('E', 129.04259), ('Q', 128.05858), ('G', 57.02146),
    ('H', 137.05891), ('I', 113.08406), ('L', 113.08406), ('K', 128.09496),
    ('M', 131.04049), ('F', 147.06841), ('P', 97.05276),  ('S', 87.03203),
    ('T', 101.04768), ('W', 186.07931), ('Y', 163.06333), ('V', 99.06841),
]
    _AA_MASS_MONO[Int(UInt8(aa)) + 1] = mass
    _AA_MASS_MONO[Int(UInt8(aa) | 0x20) + 1] = mass
end

# Average masses (Da)
for (aa, mass) in [
    ('A', 71.0788), ('R', 156.1875), ('N', 114.1038), ('D', 115.0886),
    ('C', 103.1388), ('E', 129.1155), ('Q', 128.1307), ('G', 57.0519),
    ('H', 137.1411), ('I', 113.1594), ('L', 113.1594), ('K', 128.1741),
    ('M', 131.1926), ('F', 147.1766), ('P', 97.1167),  ('S', 87.0782),
    ('T', 101.1051), ('W', 186.2132), ('Y', 163.1760), ('V', 99.1326),
]
    _AA_MASS_AVG[Int(UInt8(aa)) + 1] = mass
    _AA_MASS_AVG[Int(UInt8(aa) | 0x20) + 1] = mass
end

# ─── Hydropathicity table (Kyte-Doolittle) ────────────────────────────────────

const _HYDROPATHICITY = fill(0.0, 256)
for (aa, val) in [
    ('A', 1.8), ('R', -4.5), ('N', -3.5), ('D', -3.5), ('C', 2.5),
    ('E', -3.5), ('Q', -3.5), ('G', -0.4), ('H', -3.2), ('I', 4.5),
    ('L', 3.8), ('K', -3.9), ('M', 1.9), ('F', 2.8), ('P', -1.6),
    ('S', -0.8), ('T', -0.7), ('W', -0.9), ('Y', -1.3), ('V', 4.2),
]
    _HYDROPATHICITY[Int(UInt8(aa)) + 1] = val
    _HYDROPATHICITY[Int(UInt8(aa) | 0x20) + 1] = val
end

# ─── Amino acid code mapping (20 standard AAs → 0..19) ───────────────────────

const _AA_CODE = fill(UInt8(255), 256)
const _AA_CHARS = "ACDEFGHIKLMNPQRSTVWY"
for (index, aa) in enumerate(_AA_CHARS)
    _AA_CODE[Int(UInt8(aa)) + 1] = UInt8(index - 1)
    _AA_CODE[Int(UInt8(aa) | 0x20) + 1] = UInt8(index - 1)
end

# ─── Instability index DIWV table (20×20 matrix, indexed by _AA_CODE) ────────

const _INSTABILITY_DIWV = zeros(Float64, 20, 20)

# Fill from the Guruprasad et al. dipeptide instability weight values
const _DIWV_RAW = Dict{String,Float64}(
    "WW" => 1.0, "WC" => 1.0, "WM" => 24.68, "WH" => 24.68, "WY" => 1.0,
    "WF" => 1.0, "WQ" => 1.0, "WN" => 13.34, "WI" => 1.0, "WR" => 1.0,
    "WD" => 1.0, "WP" => 1.0, "WT" => -14.03, "WK" => 1.0, "WE" => 1.0,
    "WV" => -7.49, "WS" => 1.0, "WG" => -9.37, "WA" => -14.03, "WL" => 13.34,
    "CW" => 24.68, "CC" => 1.0, "CM" => 33.6, "CH" => 33.6, "CY" => 1.0,
    "CF" => 1.0, "CQ" => -6.54, "CN" => 1.0, "CI" => 1.0, "CR" => 1.0,
    "CD" => 20.26, "CP" => 20.26, "CT" => 33.6, "CK" => 1.0, "CE" => 1.0,
    "CV" => -6.54, "CS" => 1.0, "CG" => 1.0, "CA" => 1.0, "CL" => 20.26,
    "MW" => 1.0, "MC" => 1.0, "MM" => -1.88, "MH" => 58.28, "MY" => 24.68,
    "MF" => 1.0, "MQ" => -6.54, "MN" => 1.0, "MI" => 1.0, "MR" => -6.54,
    "MD" => 1.0, "MP" => 44.94, "MT" => -1.88, "MK" => 1.0, "ME" => 1.0,
    "MV" => 1.0, "MS" => 44.94, "MG" => 1.0, "MA" => 13.34, "ML" => 1.0,
    "HW" => -1.88, "HC" => 1.0, "HM" => 1.0, "HH" => 1.0, "HY" => 44.94,
    "HF" => -9.37, "HQ" => 1.0, "HN" => 24.68, "HI" => 44.94, "HR" => 1.0,
    "HD" => 1.0, "HP" => -1.88, "HT" => -6.54, "HK" => 24.68, "HE" => 1.0,
    "HV" => 1.0, "HS" => 1.0, "HG" => -9.37, "HA" => 1.0, "HL" => 1.0,
    "YW" => -9.37, "YC" => 1.0, "YM" => 44.94, "YH" => 13.34, "YY" => 13.34,
    "YF" => 1.0, "YQ" => 1.0, "YN" => 1.0, "YI" => 1.0, "YR" => -15.91,
    "YD" => 24.68, "YP" => 13.34, "YT" => -7.49, "YK" => 1.0, "YE" => -6.54,
    "YV" => 1.0, "YS" => 1.0, "YG" => -7.49, "YA" => 24.68, "YL" => 1.0,
    "FW" => 1.0, "FC" => 1.0, "FM" => 1.0, "FH" => 1.0, "FY" => 33.6,
    "FF" => 1.0, "FQ" => 1.0, "FN" => 1.0, "FI" => 1.0, "FR" => 1.0,
    "FD" => 13.34, "FP" => 20.26, "FT" => 1.0, "FK" => -14.03, "FE" => 1.0,
    "FV" => 1.0, "FS" => 1.0, "FG" => 1.0, "FA" => 1.0, "FL" => 1.0,
    "QW" => 1.0, "QC" => -6.54, "QM" => 1.0, "QH" => 1.0, "QY" => -6.54,
    "QF" => -6.54, "QQ" => 20.26, "QN" => 1.0, "QI" => 1.0, "QR" => 1.0,
    "QD" => 20.26, "QP" => 20.26, "QT" => 1.0, "QK" => 1.0, "QE" => 20.26,
    "QV" => -6.54, "QS" => 44.94, "QG" => 1.0, "QA" => 1.0, "QL" => 1.0,
    "NW" => -9.37, "NC" => -1.88, "NM" => 1.0, "NH" => 1.0, "NY" => 1.0,
    "NF" => -14.03, "NQ" => -6.54, "NN" => 1.0, "NI" => 44.94, "NR" => 1.0,
    "ND" => 1.0, "NP" => -1.88, "NT" => -7.49, "NK" => 24.68, "NE" => 1.0,
    "NV" => 1.0, "NS" => 1.0, "NG" => -14.03, "NA" => 1.0, "NL" => 1.0,
    "IW" => 1.0, "IC" => 1.0, "IM" => 1.0, "IH" => 13.34, "IY" => 1.0,
    "IF" => 1.0, "IQ" => 1.0, "IN" => 1.0, "II" => 1.0, "IR" => 1.0,
    "ID" => 1.0, "IP" => -1.88, "IT" => 1.0, "IK" => -7.49, "IE" => 44.94,
    "IV" => -7.49, "IS" => 1.0, "IG" => 1.0, "IA" => 1.0, "IL" => 20.26,
    "RW" => 58.28, "RC" => 1.0, "RM" => 1.0, "RH" => 20.26, "RY" => -6.54,
    "RF" => 1.0, "RQ" => 20.26, "RN" => 13.34, "RI" => 1.0, "RR" => 58.28,
    "RD" => 1.0, "RP" => 20.26, "RT" => 1.0, "RK" => 1.0, "RE" => 1.0,
    "RV" => 1.0, "RS" => 44.94, "RG" => -7.49, "RA" => 1.0, "RL" => 1.0,
    "DW" => 1.0, "DC" => 1.0, "DM" => 1.0, "DH" => 1.0, "DY" => 1.0,
    "DF" => -6.54, "DQ" => 1.0, "DN" => 1.0, "DI" => 1.0, "DR" => -6.54,
    "DD" => 1.0, "DP" => 1.0, "DT" => -14.03, "DK" => -7.49, "DE" => 1.0,
    "DV" => 1.0, "DS" => 20.26, "DG" => 1.0, "DA" => 1.0, "DL" => 1.0,
    "PW" => -1.88, "PC" => -6.54, "PM" => -6.54, "PH" => 1.0, "PY" => 1.0,
    "PF" => 20.26, "PQ" => 20.26, "PN" => 1.0, "PI" => 1.0, "PR" => -6.54,
    "PD" => -6.54, "PP" => 20.26, "PT" => 1.0, "PK" => 1.0, "PE" => 18.38,
    "PV" => 20.26, "PS" => 20.26, "PG" => 1.0, "PA" => 20.26, "PL" => 1.0,
    "TW" => -14.03, "TC" => 1.0, "TM" => 1.0, "TH" => 1.0, "TY" => 1.0,
    "TF" => 13.34, "TQ" => -6.54, "TN" => -14.03, "TI" => 1.0, "TR" => 1.0,
    "TD" => 1.0, "TP" => 1.0, "TT" => 1.0, "TK" => 1.0, "TE" => 20.26,
    "TV" => 1.0, "TS" => 1.0, "TG" => -7.49, "TA" => 1.0, "TL" => 1.0,
    "KW" => 1.0, "KC" => 1.0, "KM" => 33.6, "KH" => 1.0, "KY" => 1.0,
    "KF" => 1.0, "KQ" => 24.68, "KN" => 1.0, "KI" => -7.49, "KR" => 33.6,
    "KD" => 1.0, "KP" => -6.54, "KT" => 1.0, "KK" => 1.0, "KE" => 1.0,
    "KV" => -7.49, "KS" => 1.0, "KG" => -7.49, "KA" => 1.0, "KL" => -7.49,
    "EW" => -14.03, "EC" => 44.94, "EM" => 1.0, "EH" => -6.54, "EY" => 1.0,
    "EF" => 1.0, "EQ" => 20.26, "EN" => 1.0, "EI" => 20.26, "ER" => 1.0,
    "ED" => 20.26, "EP" => 20.26, "ET" => 1.0, "EK" => 1.0, "EE" => 33.6,
    "EV" => 1.0, "ES" => 20.26, "EG" => 1.0, "EA" => 1.0, "EL" => 1.0,
    "VW" => 1.0, "VC" => 1.0, "VM" => 1.0, "VH" => 1.0, "VY" => -6.54,
    "VF" => 1.0, "VQ" => 1.0, "VN" => 1.0, "VI" => 1.0, "VR" => 1.0,
    "VD" => -14.03, "VP" => 20.26, "VT" => -7.49, "VK" => -1.88, "VE" => 1.0,
    "VV" => 1.0, "VS" => 1.0, "VG" => -7.49, "VA" => 1.0, "VL" => 1.0,
    "SW" => 1.0, "SC" => 33.6, "SM" => 1.0, "SH" => 1.0, "SY" => 1.0,
    "SF" => 1.0, "SQ" => 20.26, "SN" => 1.0, "SI" => 1.0, "SR" => 20.26,
    "SD" => 1.0, "SP" => 44.94, "ST" => 1.0, "SK" => 1.0, "SE" => 20.26,
    "SV" => 1.0, "SS" => 20.26, "SG" => 1.0, "SA" => 1.0, "SL" => 1.0,
    "GW" => 13.34, "GC" => 1.0, "GM" => 1.0, "GH" => 1.0, "GY" => -7.49,
    "GF" => 1.0, "GQ" => 1.0, "GN" => -7.49, "GI" => -7.49, "GR" => 1.0,
    "GD" => 1.0, "GP" => 1.0, "GT" => -7.49, "GK" => -7.49, "GE" => -6.54,
    "GV" => 1.0, "GS" => 1.0, "GG" => 13.34, "GA" => -7.49, "GL" => 1.0,
    "AW" => 1.0, "AC" => 44.94, "AM" => 1.0, "AH" => -7.49, "AY" => 1.0,
    "AF" => 1.0, "AQ" => 1.0, "AN" => 1.0, "AI" => 1.0, "AR" => 1.0,
    "AD" => -7.49, "AP" => 20.26, "AT" => 1.0, "AK" => 1.0, "AE" => 1.0,
    "AV" => 1.0, "AS" => 1.0, "AG" => 1.0, "AA" => 1.0, "AL" => 1.0,
    "LW" => 24.68, "LC" => 1.0, "LM" => 1.0, "LH" => 1.0, "LY" => 1.0,
    "LF" => 1.0, "LQ" => 33.6, "LN" => 1.0, "LI" => 1.0, "LR" => 20.26,
    "LD" => 1.0, "LP" => 20.26, "LT" => 1.0, "LK" => -7.49, "LE" => 1.0,
    "LV" => 1.0, "LS" => 1.0, "LG" => 1.0, "LA" => 1.0, "LL" => 1.0,
)

for (dipeptide, value) in _DIWV_RAW
    code1 = _AA_CODE[Int(UInt8(dipeptide[1])) + 1]
    code2 = _AA_CODE[Int(UInt8(dipeptide[2])) + 1]
    code1 == 0xff && continue
    code2 == 0xff && continue
    _INSTABILITY_DIWV[Int(code1) + 1, Int(code2) + 1] = value
end

# ─── pK values for isoelectric point (IPC algorithm) ──────────────────────────

const _PK_NTERM  = 8.2     # N-terminus
const _PK_CTERM  = 3.65    # C-terminus
const _PK_D      = 3.9     # Aspartic acid
const _PK_E      = 4.07    # Glutamic acid
const _PK_C      = 8.18    # Cysteine
const _PK_Y      = 10.46   # Tyrosine
const _PK_H      = 6.04    # Histidine
const _PK_K      = 10.54   # Lysine
const _PK_R      = 12.48   # Arginine

# ══════════════════════════════════════════════════════════════════════════════
# Public API
# ══════════════════════════════════════════════════════════════════════════════

"""
    protein_mass(sequence; type="monoisotopic")

Compute the molecular weight of a protein sequence in Daltons.
The result includes the terminal water mass contribution for the peptide chain.
"""
function protein_mass(sequence::AbstractString; type::AbstractString="monoisotopic")
    mass_table = type == "monoisotopic" ? _AA_MASS_MONO : _AA_MASS_AVG
    water = type == "monoisotopic" ? 18.01524 : 18.01056
    mass = water

    @inbounds for byte in codeunits(sequence)
        residue_mass = mass_table[Int(byte) + 1]
        residue_mass == 0.0 && throw(ArgumentError("unknown amino acid '$(Char(byte))'"))
        mass += residue_mass
    end

    return mass
end

"""
    extinction_coefficient(sequence)

Compute the molar extinction coefficient at 280 nm.
Uses the Pace formula: ε = 1490×nY + 5500×nW + 125×nC.
"""
function extinction_coefficient(sequence::AbstractString)
    n_y = 0; n_w = 0; n_c = 0

    @inbounds for byte in codeunits(sequence)
        ch = byte | 0x20  # lowercase
        n_y += (ch == UInt8('y'))
        n_w += (ch == UInt8('w'))
        n_c += (ch == UInt8('c'))
    end

    return 1490 * n_y + 5500 * n_w + 125 * n_c
end

"""
    instability_index(sequence)

Compute the Guruprasad instability index for a protein sequence.
Values below 40 suggest a stable protein in vitro.
"""
function instability_index(sequence::AbstractString)
    len = ncodeunits(sequence)
    len < 2 && throw(ArgumentError("sequence must have at least 2 residues"))

    bytes = codeunits(sequence)
    diwv = _INSTABILITY_DIWV
    aa_code = _AA_CODE
    total = 0.0

    @inbounds for i in 1:(len - 1)
        code1 = aa_code[Int(bytes[i]) + 1]
        code2 = aa_code[Int(bytes[i + 1]) + 1]
        code1 == 0xff && throw(ArgumentError("unknown amino acid '$(Char(bytes[i]))'"))
        code2 == 0xff && throw(ArgumentError("unknown amino acid '$(Char(bytes[i + 1]))'"))
        total += diwv[Int(code1) + 1, Int(code2) + 1]
    end

    return 10.0 / len * total
end

"""
    isoelectric_point(sequence; precision=0.01)

Estimate the protein isoelectric point using an IPC-style bisection search.
`precision` controls the pH interval at which the search stops.
"""
function isoelectric_point(sequence::AbstractString; precision::Real=0.01)
    # Count charged residues in one pass
    n_d = 0; n_e = 0; n_c = 0; n_y = 0
    n_h = 0; n_k = 0; n_r = 0

    @inbounds for byte in codeunits(sequence)
        ch = byte | 0x20
        ch == UInt8('d') && (n_d += 1)
        ch == UInt8('e') && (n_e += 1)
        ch == UInt8('c') && (n_c += 1)
        ch == UInt8('y') && (n_y += 1)
        ch == UInt8('h') && (n_h += 1)
        ch == UInt8('k') && (n_k += 1)
        ch == UInt8('r') && (n_r += 1)
    end

    pH = 6.5
    pH_prev = 0.0
    pH_next = 14.0

    while true
        # Negative charges
        charge  = -1.0 / (1.0 + 10.0^(_PK_CTERM - pH))
        charge -= n_d / (1.0 + 10.0^(_PK_D - pH))
        charge -= n_e / (1.0 + 10.0^(_PK_E - pH))
        charge -= n_c / (1.0 + 10.0^(_PK_C - pH))
        charge -= n_y / (1.0 + 10.0^(_PK_Y - pH))
        # Positive charges
        charge += n_h / (1.0 + 10.0^(pH - _PK_H))
        charge += 1.0 / (1.0 + 10.0^(pH - _PK_NTERM))
        charge += n_k / (1.0 + 10.0^(pH - _PK_K))
        charge += n_r / (1.0 + 10.0^(pH - _PK_R))

        if charge < 0.0
            temp = pH
            pH -= (pH - pH_prev) / 2.0
            pH_next = temp
        else
            temp = pH
            pH += (pH_next - pH) / 2.0
            pH_prev = temp
        end

        (pH - pH_prev < precision) && (pH_next - pH < precision) && break
        pH >= 14.0 && break
    end

    return pH
end

"""
    gravy(sequence)

Calculate the Grand Average of Hydropathicity (GRAVY) using Kyte-Doolittle values.
"""
function gravy(sequence::AbstractString)
    len = ncodeunits(sequence)
    len == 0 && throw(ArgumentError("sequence must not be empty"))

    total = 0.0
    hydro = _HYDROPATHICITY

    @inbounds for byte in codeunits(sequence)
        total += hydro[Int(byte) + 1]
    end

    return total / len
end

"""
    aliphatic_index(sequence)

Compute the aliphatic index from the protein sequence.
"""
function aliphatic_index(sequence::AbstractString)
    len = ncodeunits(sequence)
    len == 0 && throw(ArgumentError("sequence must not be empty"))

    n_a = 0; n_v = 0; n_i = 0; n_l = 0

    @inbounds for byte in codeunits(sequence)
        ch = byte | 0x20
        n_a += (ch == UInt8('a'))
        n_v += (ch == UInt8('v'))
        n_i += (ch == UInt8('i'))
        n_l += (ch == UInt8('l'))
    end

    x_a = 100.0 * n_a / len
    x_v = 100.0 * n_v / len
    x_i = 100.0 * n_i / len
    x_l = 100.0 * n_l / len

    return x_a + 2.9 * x_v + 3.9 * (x_i + x_l)
end

"""
    protparam(sequence)

Return a bundled ProtParam-style summary of protein properties.
"""
function protparam(sequence::AbstractString)
    len = ncodeunits(sequence)
    len == 0 && throw(ArgumentError("sequence must not be empty"))

    # Single pass to count all residues and accumulate mass/hydropathicity
    mass_mono = 18.01524
    mass_avg  = 18.01056
    hydro_sum = 0.0
    n_a = 0; n_v = 0; n_i = 0; n_l = 0
    n_d = 0; n_e = 0; n_r = 0; n_k = 0
    n_c = 0; n_y = 0; n_w = 0; n_h = 0

    @inbounds for byte in codeunits(sequence)
        mass_mono += _AA_MASS_MONO[Int(byte) + 1]
        mass_avg  += _AA_MASS_AVG[Int(byte) + 1]
        hydro_sum += _HYDROPATHICITY[Int(byte) + 1]

        ch = byte | 0x20
        n_a += (ch == UInt8('a')); n_v += (ch == UInt8('v'))
        n_i += (ch == UInt8('i')); n_l += (ch == UInt8('l'))
        n_d += (ch == UInt8('d')); n_e += (ch == UInt8('e'))
        n_r += (ch == UInt8('r')); n_k += (ch == UInt8('k'))
        n_c += (ch == UInt8('c')); n_y += (ch == UInt8('y'))
        n_w += (ch == UInt8('w')); n_h += (ch == UInt8('h'))
    end

    x_a = 100.0 * n_a / len
    x_v = 100.0 * n_v / len
    x_i = 100.0 * n_i / len
    x_l = 100.0 * n_l / len

    return (
        length = len,
        molecular_weight_mono = mass_mono,
        molecular_weight_avg = mass_avg,
        negative_residues = n_d + n_e,
        positive_residues = n_r + n_k,
        extinction_coefficient = 1490 * n_y + 5500 * n_w + 125 * n_c,
        instability_index = len >= 2 ? instability_index(sequence) : NaN,
        aliphatic_index = x_a + 2.9 * x_v + 3.9 * (x_i + x_l),
        gravy = hydro_sum / len,
        isoelectric_point = isoelectric_point(sequence),
    )
end

