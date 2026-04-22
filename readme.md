# SPAR
The SPAR (Summed Product Analytic Representation) framework, is designed to serve as a general black-box interface through which kinetic energy operators (KEOs), potential energy functions (PEFs), dipole moment functions (DMFs) and possibly other operators in sum-of-products form can be communicated to a varational solver for ro-vibrational calculations of polyatomic molecules. 

SPAR is used to represent operators derived in sum-of-products form:

$`h(\mathbf{q}) = 
\sum_{{l}_1,{l}_2,\ldots,l_{3\mathcal{N}-6}}a_{{l}_1,{l}_2,\ldots,l_{3\mathcal{N}-6}} f^{(1)}_{l_1}(q_1)  f^{(2)}_{l_2}(q_2) \cdots f^{(3\mathcal{N}-6)}_{l_{3\mathcal{N}-6}}(q_{3\mathcal{N}-6})  \equiv 
\sum_{\mathbf{l}} a_{\mathbf{l}} \prod^{3\mathcal{N} - 6}_{k = 1} f^{(k)}_{l_k}(q_k).`$

Where $`h(\mathbf{q})`$ is an operator, or component of an operator, that depends on the $`3\mathcal{N}-6`$ vibrational degrees of freedom $`\mathbf{q}=(q_1, q_2,...q_{3\mathcal{N}-6})`$. Here $`\mathbf{l}`$ is a vector of indices $`\mathbf{l} = (l_1, l_2,...l_{3\mathcal{N}-6})`$ where each $`l_k`$ refers to a so-called "basic function" $`f^{(k)}_{l_k}(q_k)`$; the building blocks with which the sum-of-products form is constructed. The sum extends over all terms for which the coefficient $`a_{\mathbf{l}} \neq 0`$.

For KEOs we also include a mass-dependent form: 

$`h(\bm{q}) = 
\sum_{{l}_0,{l}_1,l_2,\ldots,} a_{{l}_0,{l}_1,l_2 \ldots,}
\frac{1}{m_{{l}_0}} f_{{l}_1}(q_1) f_{{l}_1}(q_1)   \cdots \equiv  
\sum_{\mathbf{l}} a_{\mathbf{l}}\prod^{D}_{k = 0} f^{(k)}_{{l}_k}(q_k),`$

where $`\mathbf{l}`$ is now defined such that index $`k=0`$ refers to a nuclear mass such that $`\mathbf{l} = (l_0, l_1, l_2,...l_{3\mathcal{N}-6})`$.
