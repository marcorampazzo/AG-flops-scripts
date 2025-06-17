
#library imports
import lrcalc

# my file imports
import utils

class Partition:
    """
    I made a custom class so that i can use '*' to denote the Littlewood--Richardson product
    (by the special method __mul__), to keep track of multiplicities and to have a clean print.
    Multiplicity is one by default
    """
    def __init__(self, mylist, multiplicity=1):
        self.partition = mylist
        self.multiplicity= multiplicity
    def __mul__(self, other_partition):
        raw_result = lrcalc.mult(self.partition, other_partition.partition)
        base_multiplicity = self.multiplicity * other_partition.multiplicity
        result = []
        for key, value in raw_result.items():
            result.append(Partition(list(key), value * base_multiplicity))
        return result
    def __str__(self):
        return f"{self.partition}(×{self.multiplicity})"
    def __repr__(self):
        return f"{self.partition}(×{self.multiplicity})"

class HomogeneousIrreducible:
    """
    this class represents homogeneous irrediucible vector bundles on G(k, n).
    the first partition is a list representing a Schur power of U*, while
    the second represents a Schur power of Q
    Multiplicity is one by default
    """
    def __init__(self, k, n, first_partition, second_partition, multiplicity=1):
        self.k = k
        self.n = n
        """
        beware that self.first_partition is a Partition object. To access the partition as a list,
        call self.first_partition.partition. Same for second_partition
        """
        self.first_partition = Partition(first_partition + [0] * (k - len(first_partition)), 1)
        self.second_partition = Partition(second_partition + [0] * (n - k - len(second_partition)), 1)
        self.multiplicity= multiplicity
    def __mul__(self, other_double_partition):
        """
        we compute the product of two HomogeneousVectorBundle's with '*'. Littlewood--Richardson is applied
        on both partitions.
        """
        result = []
        result_first_partitions = self.first_partition * other_double_partition.first_partition
        result_second_partitions = self.second_partition * other_double_partition.second_partition
        for item in result_first_partitions:
            for other_item in result_second_partitions:
                if len(item.partition) <= self.k and len(other_item.partition) <= self.n-self.k:
                    total_multiplicity = self.multiplicity * other_double_partition.multiplicity * item.multiplicity * other_item.multiplicity
                    result.append(HomogeneousIrreducible(self.k, self.n, item.partition, other_item.partition, total_multiplicity))
        return HomogeneousDirectSum(result)
    def __str__(self):
        return f"{self.first_partition.partition}|{self.second_partition.partition}(×{self.multiplicity})"
    def __repr__(self):
        return f"{self.first_partition.partition}|{self.second_partition.partition}(×{self.multiplicity})"
    def rank(self):
        return utils.weyl_dim(self.first_partition.partition) * utils.weyl_dim(self.second_partition.partition)
    def dual(self):
        # reverse and negate
        part1 = self.first_partition.partition
        part2 = self.second_partition.partition
        dual_first = [-x for x in reversed(part1)]
        dual_second = [-x for x in reversed(part2)]
        # find the global twist
        max_entry = -min(dual_first + dual_second)
        t = max_entry
        # shift by k
        twisted_first = [x + t for x in dual_first]
        twisted_second = [x + t for x in dual_second]
        # return a new HomogeneousIrreducible (with same multiplicity)
        return HomogeneousIrreducible(self.k, self.n, twisted_first, twisted_second, self.multiplicity)

    def cohomology(self):
        """
        compute cohomology with Borel--Weyl--Bott.
        If rep is true, the output has the form {"representation": [,,,], "dimension": k, "degree": p},
        otherwise  {"representation": [,,,], "degree": p}
        """
        part1 = self.first_partition.partition
        part2 = self.second_partition.partition
        # the following two lines are reduntant I think...
        padded_first = part1 + [0] * (self.k - len(part1))
        padded_second = part2 + [0] * (self.n - self.k - len(part2))
        total_partition = padded_first + padded_second
        # adding rho
        total_partition = utils.add_rho(total_partition)
        # checking repetitions.
        # if there are none:
        if not utils.has_repeats(total_partition):
            degree = 0
            while True:
                if utils.is_decreasing(total_partition):
                    total_partition = utils.subtract_rho(total_partition)
                    return {
                        "representaiton": utils.normalize_partition(total_partition),
                        "dimension": utils.weyl_dim(total_partition),
                        "degree": degree
                    }
                else:
                    utils.swap_first_increase(total_partition)
                    degree = degree + 1
        # acyclic case
        else:
            return "no cohomology"


class HomogeneousDirectSum(list):
    """
    this class is just a wrapper for HomogeneousIrreducible. I made it so that
    the tensor product is distributive wrt direct sum (as it is supposed to be)
    """
    def __init__(self, summands=None):
        """
        summands: an iterable of HomogeneousIrreducible (or other DirectSums).
        """
        # if they passed nothing, treat as empty list
        super().__init__(summands or [])
    def __mul__(self, other):
        """
        Distribute self * other over every summand.
        other may be a single HI or another direct sum.
        """
        out = []
        # if other is a “vector”:
        if isinstance(other, HomogeneousDirectSum):
            for A in self:
                for B in other:
                    # A*B returns another direct sum
                    prod = A * B
                    out.extend(prod if isinstance(prod, list) else [prod])
        else:
            # other is a single HI
            for A in self:
                prod = A * other
                out.extend(prod if isinstance(prod, list) else [prod])
        return HomogeneousDirectSum(out)

    def dual(self):
        """
        Return the direct sum of the duals of each summand.
        """
        return HomogeneousDirectSum([X.dual() for X in self])

    def __rmul__(self, other):
        """
        Called when the left operand doesn't know how to multiply
        by a direct sum.  Typically other is a single HI.
        """
        # other * self  = distribute other over our list
        out = []
        for A in self:
            prod = other * A
            out.extend(prod if isinstance(prod, list) else [prod])
        return HomogeneousDirectSum(out)

    def cohomology(self):
        """
        cohomology of a direct sum is the direct sum of the cohomologies
        """
        return [X.cohomology() for X in self]

