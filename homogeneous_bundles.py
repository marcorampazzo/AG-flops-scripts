
#library imports
import lrcalc

# my file imports
import utils

class Partition:
    """
    We made a custom class so that we can use '*' to denote the Littlewood--Richardson product
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
    Twist is a feature of this specific problem. If twist is nonzero, this class describes a pullback of the bundle
    to the projectivization of Q*+O(-1) on G(k, n), twisted by O_rel (twist).
    """
    def __init__(self, k, n, first_partition, second_partition, multiplicity=1, twist=0, shift=0):
        self.k = k
        self.n = n
        """
        beware that self.first_partition is a Partition object. To access the partition as a list,
        call self.first_partition.partition. Same for second_partition
        """
        self.first_partition = Partition(first_partition + [0] * (k - len(first_partition)), 1)
        self.second_partition = Partition(second_partition + [0] * (n - k - len(second_partition)), 1)
        self.multiplicity= multiplicity
        self.twist = twist
        self.shift = shift
    def __mul__(self, other_irreducible_bundle):
        """
        we compute the product of two HomogeneousVectorBundle's with '*'. Littlewood--Richardson is applied
        on both partitions.
        """
        result = []
        if isinstance(other_irreducible_bundle, HomogeneousIrreducible):
            result_first_partitions = self.first_partition * other_irreducible_bundle.first_partition
            result_second_partitions = self.second_partition * other_irreducible_bundle.second_partition
            for item in result_first_partitions:
                for other_item in result_second_partitions:
                    if len(item.partition) <= self.k and len(other_item.partition) <= self.n-self.k:
                        total_multiplicity = self.multiplicity * other_irreducible_bundle.multiplicity * item.multiplicity * other_item.multiplicity
                        total_twist = self.twist + other_irreducible_bundle.twist
                        total_shift = self.shift + other_irreducible_bundle.shift
                        result.append(HomogeneousIrreducible(self.k, self.n, item.partition, other_item.partition, total_multiplicity, total_twist, total_shift))
        elif isinstance(other_irreducible_bundle, HomogeneousDirectSum):
            result = []
            for summand in other_irreducible_bundle:
                result.append(self * summand)
            result = utils.flatten(result)
        return HomogeneousDirectSum(result)
    def __str__(self):
        return (
            f"{self.first_partition.partition}|{self.second_partition.partition}"
            f" × O_rel({self.twist})"
            f"[{self.shift}]"
            f"(×{self.multiplicity})"
    )
    def __repr__(self):
        return (
            f"{self.first_partition.partition}|{self.second_partition.partition}"
            f" × O_rel({self.twist})"
            f"[{self.shift}]"
            f"(× {self.multiplicity})"
    )
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
        return HomogeneousIrreducible(self.k, self.n, twisted_first, twisted_second, self.multiplicity, -self.twist)

    def pushforward(self):
        if self.twist == 0:
            return self
        elif self.twist > 0:
            list_of_bundles = [] # the summands appearing  in the pushforward of O(twist)
            for l in range(self.twist+1):
                list_of_bundles.append(HomogeneousIrreducible(self.k, self.n, [self.twist for i in range(self.k)], [l for i in range(self.n - self.k - 1)], self.shift))
            pushforward_of_lb = HomogeneousDirectSum(list_of_bundles)
                # print(f"pushforward of the line bundle: {pushforward_of_lb}")
            self.twist = 0
            return pushforward_of_lb * self
        elif self.twist > -5:
            return HomogeneousIrreducible(self.k, self.n, [], [], 0, 0, 0) # the zero bundle is O with multiplicity 0 (i know it's horrible...)
        else: # we use Serre duality
            canonical = HomogeneousIrreducible(self.k, self.n, [], [4 for i in range(self.n - self.k)], 1, -5)
            new_bundle = self.dual() * canonical
            if isinstance(new_bundle, HomogeneousIrreducible): 
                new_bundle.shift = new_bundle.shift + 4
            elif isinstance(new_bundle, HomogeneousDirectSum):
                for item in new_bundle:
                    item.shift = item.shift + 4
            return new_bundle.pushforward()

    def bott(self):
        """
        compute cohomology with Borel--Weyl--Bott.
        If rep is true, the output has the form {"representation": [,,,], "dimension": k, "degree": p},
        otherwise  {"representation": [,,,], "degree": p}
        """
        if self.multiplicity == 0:
            return 'acyclic'
        part1 = self.first_partition.partition
        part2 = self.second_partition.partition
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
                        "representation": utils.normalize_partition(total_partition),
                        "dimension": utils.weyl_dim(total_partition),
                        "degree": degree + self.shift,
                        "multiplicity": self.multiplicity
                    }
                else:
                    utils.swap_first_increase(total_partition)
                    degree = degree + 1
        # acyclic case
        else:
            return "acyclic"
        
    def cohomology(self):
        if self.multiplicity == 0:
            return 'acyclic'
        return self.pushforward().bott()

    def euler(self):
        coh = self.cohomology()
        # print(coh)
        # if coh[0] == 'acyclic':
        #     return 0

        if coh == 'acyclic':
            return 0
        
        if isinstance(coh, dict):
            # print(coh)
            return ((-1)**(coh['degree'])) * coh['dimension'] * coh['multiplicity']
        else:
            out = 0
            for item in coh:
                # print(item)
                if item != 'acyclic':
                    out = out + ((-1)**(item['degree'])) * item['dimension'] * coh['multiplicity']
            return out


class HomogeneousDirectSum(list):
    """
    this class is just a wrapper for HomogeneousIrreducible. We made it so that
    the tensor product is distributive wrt direct sum
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
    
    def pushforward(self):
        """
        pushforward of a direct sum is the direct sum of the pushforwards
        """
        out = [X.pushforward() for X in self]
        out = utils.flatten(out)
        return HomogeneousDirectSum(out)
    
    def bott(self):
        """
        cohomology of a direct sum is the direct sum of the cohomologies
        """
        out = [X.bott() for X in self]
        return utils.flatten(out)
    
    def cohomology(self):
        """
        cohomology of a direct sum is the direct sum of the cohomologies
        """
        out = [X.cohomology() for X in self]
        return utils.flatten(out)

    def euler(self):
        """
        Euler characteristic of a direct sum is the direct sum of the Euler charactieristics
        """

        out = 0

        for X in self:
            out = out + X.euler() 
        
        return out
