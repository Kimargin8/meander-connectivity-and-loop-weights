"""Count all windows by interval decomposition and reflection, without listing them."""
from collections import defaultdict
from functools import cache


def oriented_with_split(total, cut, depth):
    left, right = max(0, cut-depth), min(total, cut+depth)

    @cache
    def interval(start, stop):
        if start == stop:
            return 1
        result = 0
        for partner in range(start+1, stop, 2):
            if not (left <= start < right or left <= partner < right):
                continue
            classes = 2 if start < cut <= partner else 1
            result += classes*interval(start+1, partner)*interval(partner+1, stop)
        return result

    return interval(0, total)


def fixed_at_midpoint(half, depth):
    # The first half of a symmetric matching is a Dyck prefix. Every far
    # opening must be closed on this same side. An unmatched end meets its
    # mirror and has two choices of class. Input order is outside-in.
    far = max(0, half-depth)
    paths = {(far, far): 1}
    for _ in range(far, half):
        updated = defaultdict(int)
        for (height, far_left), count in paths.items():
            updated[height+1, far_left] += count
            if height:
                updated[height-1, min(height-1, far_left)] += count
        paths = updated
    return sum(2**h * count for (h, far_left), count in paths.items() if not far_left)


def count_windows(depth):
    oriented = fixed = 0
    for total in range(0, 4*depth+1, 2):
        oriented += sum(oriented_with_split(total, cut, depth) for cut in range(total+1))
        fixed += fixed_at_midpoint(total//2, depth)
    assert (oriented+fixed) % 2 == 0
    return dict(oriented=oriented, fixed=fixed, orbits=(oriented+fixed)//2)


if __name__ == "__main__":
    for d in range(1, 9):
        print(d, count_windows(d))
