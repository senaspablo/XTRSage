# XTRSage
This is a humble implementation of XTR (Efficient and Compact Subgroup Trace Representation) using SageMath.

Broadly speaking, the XTR system is based on the discrete logarithm in a cyclic group, but instead using (field) traces of powers of one of the generators. It allows us to use the arithmetic of ![$GF(p^2)$](https://render.githubusercontent.com/render/math?math=%24GF(p%5E2)%24), which allows a very efficient implementation, while working with elements belonging to ![$GF(p^6)$](https://render.githubusercontent.com/render/math?math=%24GF(p%5E6)%24).

More information: [An Overview of the XTR Public Key System (2000) ](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.104.2847)
