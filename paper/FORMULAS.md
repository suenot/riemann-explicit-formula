# Formulas and numerical conventions

All logarithms are natural. The cached positive ordinates $\gamma_k$ are
ordered by increasing height. A budget $N$ includes $\rho_k=\tfrac12+i\gamma_k$
and its conjugate for each $k\le N$. E1–E5 use binary64 arithmetic for the
sampled sums. E1's wave-extremum check uses 38 decimal places and E6 uses 40
with higher precision for Gram coefficients and its independent check.

## E1. Waves

The arithmetic functions use $\Lambda(p^a)=\log p$ for a prime $p$ and
integer $a\ge1$, and $\Lambda(n)=0$ otherwise. The staircase is
$\psi(x)=\sum_{n\le x}\Lambda(n)$. Its midpoint $\psi_0(x)$ gives the term at
$n=x$ half weight. The independent right-continuous check uses
$\log\operatorname{lcm}(1,\ldots,\lfloor x\rfloor)$, averaged with the
preceding integer's value when $x$ is integral.

For $x>1$, the finite explicit formula with weights $w_k$ is

$$
\psi_{N,w}(x)=x-\log(2\pi)-\tfrac12\log(1-x^{-2})
 -\sum_{k=1}^N w_k A_k(x)\cos(\gamma_k\log x-\arg\rho_k),
\qquad A_k(x)=\frac{2\sqrt{x}}{\sqrt{1/4+\gamma_k^2}}.
$$

The displayed wave equals $2\Re(x^{\rho_k}/\rho_k)$. Its period in $\log x$
is $2\pi/\gamma_k$. For the independent complex-power check, find
successive stationary points of
$2\Re(e^{\rho_k u}/\rho_k)e^{-u/2}$ starting at the predicted maximum,
minimum and next maximum. The distance between the maxima measures the
period. Half the peak-to-trough difference, multiplied by $\sqrt{20}$,
measures the envelope at $x=20$.

E1 takes $w_k=1$. The worked phase example uses $\gamma=14.1347$ and
reduces $\gamma\log x-\arg\rho$ modulo $2\pi$. For each set of points,
$L^1$ is $\sum_i|\psi_{N,1}(x_i)-\psi_0(x_i)|$, MAE is that sum
divided by the number of points, and maximum error is the largest summand.
For $N>0$, the reference truncation scale is
$x\log^2(xT)/T+\log x\min(1,x/(Td))$, with $T=\gamma_N+10^{-12}$ and
$d$ the distance to the closest prime power. The scale has an unspecified
big-O constant. [Kedlaya's notes](https://kskedlaya.org/ant/chap-von-mangoldt.html).

## E2. Gibbs phenomenon and smoothing

For $r=\gamma_k/\gamma_N$, the weights are raw $1$, Fejér $1-r$, Gaussian
$e^{-r^2}$ and Hann $(1+\cos(\pi r))/2$. The Gaussian is truncated at the
same last ordinate as all the other windows. All weights are empty at $N=0$.

At a prime $p$, with $q$ the next prime power, the right-side interval in
$u=\log x$ is $[\log p,\min((\log p+\log q)/2,\log p+2\pi/\gamma_N)]$.
The maximum is taken over both endpoints and over positive-to-negative
crossings of

$$
\frac{d}{du}\psi_{N,w}(e^u)
=e^u-\frac{1}{e^{2u}-1}-2e^{u/2}\sum_{k=1}^N w_k\cos(\gamma_k u).
$$

The derivative is scanned with 48 intervals and rechecked with 96.
The signed excess is the maximum minus the right-continuous plateau
$\psi(p)$. The overshoot percent is 100 times its positive part divided
by $\log p$. The classical Fourier benchmark is
$100(\operatorname{Si}(\pi)/\pi-1/2)$ with
$\operatorname{Si}(z)=\int_0^z\sin(v)/v\,dv$. Tapers are scored by MAE
on $x=1.5,1.55,\ldots,50$ and worst overshoot at $p=7,11,13$. Selection
uses the smallest tested budget with MAE $\le0.25$ and overshoot $\le2\%$,
with lowest MAE breaking ties.

## E3. Prime-power spectrum

Let $a=\log 2$, $L=\log(X/2)$ and
$W(u)=(1-\cos(2\pi(u-a)/L))/2$ on $[a,a+L]$, zero elsewhere. The
signed windowed cosine sum and the amplitude used for peak search are

$$
P_X(t)=\sum_{n\le X}\frac{\Lambda(n)}{\sqrt n}W(\log n)\cos(t\log n),
\qquad
S_X(t)=\left|P_X(t)-\int_a^{a+L}e^{u/2}W(u)\cos(tu)\,du\right|.
$$

The integral is the same projection of $d(e^u)$ after multiplication by
$e^{-u/2}$. Put $z=\tfrac12+it$, $\omega=2\pi/L$ and
$I(z)=(e^{zL}-1)/z$. The integral is the real part of
$e^{za}[I(z)/2-I(z+i\omega)/4-I(z-i\omega)/4]$.

Use $t=0,0.02,\ldots,60$. Select the ten largest local maxima of $S_X$
in $[10,50]$ by amplitude, then fit a parabola to each grid maximum and
its two neighbors. Only after selection, match each refined peak to the
nearest reference ordinate if its separation is at most half the
effective spacing $2\pi/L$. Keep at most one match per zero. The first
peak's full width at half height uses linear interpolation on each
side of the largest grid peak in $[10,18]$. The spacing scales are
$2\pi/\log X$ and $2\pi/L$. The measured width is a property of this
real cosine projection and its phase. Resolving nearby lines is a
separate question.

For the parabolic fit in E3 and E4, let $y_-,y_0,y_+$ be consecutive
samples spaced by $h$. With $d=y_--2y_0+y_+$,
the vertex is displaced by $h\delta$, where $\delta=(y_--y_+)/(2d)$,
and its height is $y_0-(y_--y_+)\delta/4$.

## E4. Prime powers from zeros

Compute $F(t)=\sum_{k=1}^{1000}W_k\cos(\gamma_kt)$ for
$t=0.3,0.302,\ldots,3.5$, with $G=\gamma_{1000}$ and
$W_k=(1+\cos(\pi\gamma_k/G))/2$. For each integer $2\le n\le31$,
find the lowest interpolated minimum within $\pi/G$ of $\log n$.
Record both $F(\log n)$ and the minimum when present. Location error is
$|t_{\min}-\log n|$.

For depth, evaluate $F(\log n+d)$ at
$d\in\{-0.06,-0.05,-0.04,0.04,0.05,0.06\}$. The depth is the fourth-smallest
of these six values minus $F(\log n)$.

Landau's fixed-$n$ main term gives
$\sum_{0<\gamma\le T}\cos(\gamma\log n)\sim
 -T\Lambda(n)/(2\pi\sqrt n)$ when $\Lambda(n)>0$ and the zeros in the
sum have real part $\tfrac12$. More generally the displayed expression
is the linear-in-$T$ term, which vanishes at non-prime-powers.
Partial summation with the Hann-type weight replaces $T$ by
$\int_0^G(1+\cos(\pi y/G))/2\,dy=G/2$. The predicted positive depth at
$p^a$ is therefore $\frac{G}{4\pi}\frac{\log p}{\sqrt n}$.
The code records depth times $\sqrt n$ and divides it by $G/(4\pi)$
to obtain the normalized depth $D(n)$. At each prime power, the relative
error is $|D(n)/\log p-1|$. A pronounced dip has $D(n)>0.3$.
[Landau (1912)](https://link.springer.com/article/10.1007/BF01456808).

## E5. Bohr coefficients

The absolutely convergent Dirichlet series on $\Re s=2$ gives the
finite-time average by termwise integration,

$$
C_T(n)=\frac{1}{T}\int_0^T\zeta(2+it)n^{it}\,dt
 =\sum_{m=1}^\infty m^{-2}e^{i\theta_m}\frac{\sin\theta_m}{\theta_m},
\qquad \theta_m=\tfrac T2\log(n/m).
$$

At $m=n$ the quotient is 1. The calculation truncates at $M=100000$.
For $n\le12$ the omitted tail is bounded in modulus by
$2/(TM\log(M/12))$, since $\sum_{m>M}m^{-2}\le1/M$.
The limit is $n^{-2}$ and the multiplicative check uses
$6^{-2}=2^{-2}3^{-2}$. Finite-time errors are the complex moduli
$|C_T(n)-n^{-2}|$ and $|C_T(6)-C_T(2)C_T(3)|$.

## E6. Prime counting

Use $\pi(x)=\sum_{p\le x}1$ and $J(x)=\sum_{p^a\le x}1/a$.
The midpoint versions $\pi_0,J_0$ halve the contribution at equality.
For either matching convention,

$$
\pi(x)=\sum_{m\le\log_2 x}\frac{\mu(m)}{m}J(x^{1/m}),
\qquad
\pi_0(x)=\sum_{m\le\log_2 x}\frac{\mu(m)}{m}J_0(x^{1/m}).
$$

Here $\mu(1)=1$, $\mu(m)=0$ if $m$ has a squared prime factor and
$\mu(m)=(-1)^k$ for a product of $k$ distinct primes. Exact comparisons
at roots use $p^{am}\le x$ and $p^{am}=x$ with integer powers. $J$
and its arithmetic inversion use rational arithmetic.

The logarithmic integral is $\operatorname{Li}(x)=\operatorname{Ei}(\log x)$,
using the principal-value integral from zero and no offset at $x=2$.
For $x>1$ the real Riemann $R(x)$ is Gram's function evaluated at $\log x$.
The entire function of the log argument is

$$
G_R(z)=1+\sum_{k\ge1}\frac{z^k}{k\,k!\,\zeta(k+1)}.
$$

The finite diagnostic is
$H_N(x)=R(x)-2\sum_{k\le N}\Re G_R(\rho_k\log x)$, with no reduction
of the imaginary part of $\rho_k\log x$. For the complex evaluation,
set $b_k=1/\zeta(k+1)-\sum_{m\le B}\mu(m)/m^{k+1}$ and use

$$
G_R(z)=1+
\sum_{m\le B}\frac{\mu(m)}m
 \left(\operatorname{Ei}(z/m)-\gamma_{\mathrm E}-\log(z/m)\right)
+\sum_{k\ge1}\frac{b_k z^k}{k\,k!}.
$$

The split uses $B=100$ (rechecked at 120), $k\le100$ and
240-decimal coefficient evaluation. The arguments have positive real
part, where $\operatorname{Ei}(w)-\gamma_{\mathrm E}-\log w$ equals the
entire series $\sum_{k\ge1}w^k/(k\,k!)$ using the principal log.
For the first omitted index $K=101$ and $r=|z|/B<K+1$, the remaining
series is bounded by $r^K/[K^2 K!(1-r/(K+1))]$. This follows from
$|b_k|\le1/(kB^k)$ and a geometric bound on consecutive terms.
Grobner proves divergence with ordinary height ordering.
[Grobner (2026), v2](https://arxiv.org/abs/2609.02713v2).

The convergent reconstruction forms the explicit-formula truncation
for $y>1$,

$$
J_N(y)=\operatorname{Li}(y)-
2\sum_{k\le N}\Re\operatorname{Ei}(\rho_k\log y)-\log 2+
\int_y^\infty\frac{dt}{t(t^2-1)\log t},
$$

and then uses $M_N(x)=\sum_{m\le\log_2 x}\mu(m)J_N(x^{1/m})/m$.
With the full height-ordered zeros and their actual real parts, the
explicit formula for $J$ gives the midpoint $J_0$ and hence
$M_N(x)\to\pi_0(x)$ at fixed $x$. The code uses the finite cached
critical-line ordinates. Signed errors are estimate minus the matching
arithmetic count.
