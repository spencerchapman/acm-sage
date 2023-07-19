from collections import Counter


class ArithmeticalCongruenceMonoid:

    def __init__(self, a, b):
        if ((a ^ 2) % b) != (a % b):
            raise ValueError("a and a^2 must be equivalent mod b")

        self.b = b
        self.a = a % b

        self.__factorizations = {1: [[]]}
        self.__irreducibles = {1: False}
        self.__lengthsets = {1: [0]}

    def Factorizations(self, num):
        if num in self.__factorizations:
            return self.__factorizations[num]
        divs = [i for i in divisors(num) if i in self and num / i in self]
        return self.__ACMFactor(num, divs)

    def __ACMFactor(self, num, divs):
        for d in divs:
            if d in self.__factorizations:
                continue
            self.__factorizations[d] = []
            for s in divs:
                if s >= d:
                    break
                if self.__irreducibles[s] and d / s in self:
                    for f in self.__factorizations[(d / s)]:
                        if s >= f[-1]:
                            self.__factorizations[d].append(f + [s])
            if len(self.__factorizations[d]) == 0:
                self.__factorizations[d] = [[d]]
                self.__irreducibles[d] = True
            else:
                self.__irreducibles[d] = False
        return self.__factorizations[num]

    def NumberOfFactorizations(self, num):
        return len(self.Factorizations(num))

    def FactorizationsUpToElement(self, nmax):
        # Finds all arithmetic factorizations up to a certain element. Useful to run before long calculations
        for i in range(self.a, nmax + 1, self.b):
            self.Factorizations(i)

    def IsIrreducible(self, num):
        # checks if an element is irreducible.
        if num in self.__irreducibles:
            return self.__irreducibles[num]

        if num == 1:
            return False

        if num in self.__factorizations:
            return max(self.LengthSet(num)) == 1

        # removes all numbers not divisors in the monoid
        divs = [i for i in divisors(num) if i in self and num / i in self]

        self.__irreducibles[num] = (len(divs) <= 2)
        return self.__irreducibles[num]

    def LengthSet(self, num, t=None):
        if t is None or t == 1:
            return self.__ACMLengthSet(num)
        if t == 0:
            return set(len(set(F)) for F in self.Factorizations(num))
        if t == oo:
            return set(Counter(F).most_common()[0][1] for F in self.Factorizations(num))
        raise ValueError("t = ", t, "not implemented")

    def __ACMLengthSet(self, num):
        if num in self.__lengthsets:
            return self.__lengthsets[num]
        divs = [i for i in divisors(num) if i in self and num / i in self]
        for d in divs:
            if d in self.__lengthsets:
                continue
            self.__lengthsets[d] = []
            for s in divs:
                if s >= d:
                    break
                if self.__irreducibles[s] and d / s in self:
                    for l in self.__lengthsets[(d / s)]:
                        self.__lengthsets[d].append(l + 1)
            if len(self.__lengthsets[d]) == 0:
                self.__lengthsets[d] = [1]
                self.__irreducibles[d] = True
            else:
                self.__irreducibles[d] = False
            self.__lengthsets[d] = list(set(self.__lengthsets[d]))
        return self.__lengthsets[num]

    def MaxFactorizationLength(self, num, t=None):
        # finds the max factorization length of a given element
        return max(self.LengthSet(num, t))

    def MinFactorizationLength(self, num, t=None):
        # finds the min factorization length of a given element
        return min(self.LengthSet(num, t))

    def Elasticity(self, num, t=None):
        # finds elasticity of an element: max factorization length / min...
        if num == 1:
            return 1
        return float(self.MaxFactorizationLength(num, t)) / float(self.MinFactorizationLength(num, t))

    @staticmethod
    def MonoidsForGivenB(b):
        return [a for a in range(b) if (a ^ 2 % b) == a]

    def __eq__(self, other):
        return (self.a, self.b) == (other.a, other.b)

    def __ne__(self, other):
        return (self.a, self.b) != (other.a, other.b)

    def __contains__(self, other):
        return int(other) == other and (other == 1 or other % self.b == self.a)

    def __repr__(self):
        return "Arithmetical Congruence Monoid " + str((self.a, self.b))

    def __str__(self):
        return self.SaveToString()
