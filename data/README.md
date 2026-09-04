# Data README

These files have been produced with XATOM as follows, for example:

```
# bound states
xatom -s Ne -bound | sed -e 's/^#[[:space:]]*r(a.u.)/r(a.u.)/' > Ne.txt

# unbound states
for e in {950..1050..1}
do
    xatom -s Ne -continuum $e | sed -e 's/^#[[:space:]]*r(a.u.)/r(a.u.)/' > Ne_${e}.txt
done

# photo-ionization cross section
xatom -s Ne -pcs -PE 0-1500 -dE 1 | sed -e 's/^#[[:space:]]*P.E.(eV)/P.E.(eV)/' > pcs_Ne.txt
```

The sed expression removes the comment from the header line, so that it can be read too, to identify the orbitals.

