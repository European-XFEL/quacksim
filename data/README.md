# Data README


These files have been produced with XATOM as follows, for example, for Ne:

```
# photo-ionization cross section
xatom -s Ne -pcs -PE 0-3000 -dE 1 | sed -e 's/^#[[:space:]]*P.E.(eV)/P.E.(eV)/' > pcs_Ne.txt

# with transition-dipole moment
xatom -s Ne -pcs -PE 0-3000 -dE 1 -v | sed -e 's/^#[[:space:]]*P.E.(eV)/P.E.(eV)/' > pcs_Ne_v.txt
```

The sed expression removes the comment from the header line, so that it can be read too, to identify the orbitals.

After the files are produced, they can be read and converted into an HDF5 file for fast access using the notebook in `notebooks/XATOM_cross_sections_to_H5.ipynb`.

