# Classical & quantum security levels for generalized XMSS + Hashing-into-Hypercubes (with rejection sampling).
# Based on [DKKW25c](https://eprint.iacr.org/2025/055.pdf) Sec. 6 and [HKKTW26](https://eprint.iacr.org/2026/016) Corollary 1

import math

v = 46  # number of Winternitz chunks
w = 3  # chunk bit width; each chunk in {0, …, 2^w − 1}
pub_param = 5 * 31  # public-parameter length log₂|P| in bits
randomness = 7 * 31  # randomness / seed length log₂|R| in bits
digest = 8 * 31  # chain / tree hash digest log₂|H| in bits
lifetime = 1 << 32  # number of signing epochs
K = 100_000  # max encoding retries (target-sum Winternitz)
prime = 2**31 - 2**24 + 1  # KoalaBear prime; used for abort correction from [HKKTW26]

# ---------------------------------------------------------------------------

assert 24 % w == 0, "w must divide 24 for the rejection sampling to work nicely with koalabear"
z = 24 // w # use for rejection sampling; see [HKKTW26] Sec. 6.1

# Bit sizes of the parameter spaces
bits_H = digest
bits_P = pub_param
bits_R = randomness
bits_msg = v * w

# Abort correction (Section 6.1 of [HKKTW26])
wz = (1 << w)**z
Q = prime // wz
ell = math.ceil(v / z)

non_abort_total = ((Q * wz) / prime) ** ell
abort_correction_bits = -math.log2(non_abort_total)

bits_msg_eff = bits_msg + abort_correction_bits

# Derived quantities
log5 = math.log2(5)
log12 = math.log2(12)
log3 = math.log2(3)
logL = math.log2(lifetime)
logv = math.log2(v)
logK = math.log2(K)
logqs = math.log2(lifetime)

# Classical: each constraint gives an upper bound on k_C
k_C = min(
    bits_H - log5 - 2 * w - logL - logv,  # digest  (SM-TCR/UD/PRE)
    bits_P - log5 - 3,  # pub_param (SM-TCR)
    bits_msg_eff - log5 - 1,  # msg_hash (rTCR coll)
    bits_R - log5 - logqs - logK - 1,  # randomness (rTCR reprog)
)

# Quantum: each constraint gives an upper bound on k_Q
k_Q = min(
    bits_H / 2 - log5 - 2 * w - logL - logv - log12,  # digest
    (bits_P - 5) / 2 - log5 - 2,  # pub_param
    (bits_msg_eff - 3) / 2 - log5 - 1,  # msg_hash
    (bits_R - logqs) / 2 - log5 - log3 - logK,  # randomness
)

print(f"classical = {k_C:.2f}")
print(f"quantum   = {k_Q:.2f}")
