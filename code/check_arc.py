"""Certify the contour |z| = 3/2 for the profile theorem on the phase arc pi/7 <= theta <= pi."""
from math import factorial

import whittaker_certificate as certificate


class CertificateFailure(ArithmeticError):
    """Raised when a certified inequality fails; never disabled by python -O."""


def require(condition, message='certificate inequality failed'):
    if not condition:
        raise CertificateFailure(message)



def main():
    radius = certificate.iv.mpf(3) / 2
    degree = certificate.DEGREE
    # Every successive tail term has ratio at most (3/2)/(degree + 2).
    certificate.TAIL = 2 * radius ** (degree + 1) / factorial(degree + 1)
    power_index = len(certificate.POWERS)
    certificate.POWERS.append([radius**k / factorial(k) for k in range(degree + 1)])
    certificate.box_error.cache_clear()

    pending = [(0, 0, 0)]
    excluded = skipped = maximum_depth = 0
    while pending:
        level, i, j = pending.pop()
        denominator = 2 ** (level + 1)
        # The theta interval is [2i, 2i+2] * pi / denominator.
        if 7 * (2*i + 2) < denominator:
            skipped += 1
            continue
        if certificate.excludes(level, i, j, power_index):
            excluded += 1
            maximum_depth = max(maximum_depth, level)
            continue
        if level >= 22:
            raise ArithmeticError("Failed to certify the phase arc at the prescribed depth.")
        pending.extend((level + 1, 2*i + a, 2*j + b) for a in (0, 1) for b in (0, 1))
    print("No zero on |z|=3/2 for pi/7 <= theta <= pi.")
    print("Excluded boxes:", excluded, "boxes below lower angle:", skipped,
          "maximum depth:", maximum_depth)


if __name__ == "__main__":
    main()
