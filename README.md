# acm-sage
A fork of [acm-sage]('https://github.com/coneill-math/acm-sage').
## Additional Features
* Dynamic computation of factorizations
* Factorization Length computation under various norms.
## Usage
Download `ArithmeticalCongruenceMonoid.sage` and load the filepath:
```python
load('path/to/file/ArithmeticalCongruenceMonoid.sage')
```
A demonstration of the `ArithmeticalCongruenceMonoid` class methods:
```python
load('path/to/file/ArithmeticalCongruenceMonoid.sage')
hilbert = ArithmeticalCongruenceMonoid(1,4)
print hilbert.Factorizations(441)

print hilbert.LengthSet(441,1)
print hilbert.LengthSet(441,oo)
```