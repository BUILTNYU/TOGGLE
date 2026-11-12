# Optimal mode switching under an additive Ornstein–Uhlenbeck process

This is the entry–exit ("market switching") model of Dixit (1989), solved for demand that follows an additive Ornstein–Uhlenbeck process. It corrects the derivation in Guo, Chow & Schonfeld (2018), which used the proportional-volatility diffusion $\sigma Q\,dW$ in the valuation but simulated additive noise (see the last section).

## 1. Model

Demand follows

$$
dQ_t=\mu\,(m-Q_t)\,dt+\sigma\,dW_t,\qquad \mu>0,\ \sigma>0 .
$$

Given $Q_0=Q$, $Q_s\sim N\!\big(m+(Q-m)e^{-\mu s},\ \tfrac{\sigma^2}{2\mu}(1-e^{-2\mu s})\big)$, and the long-run distribution is $N(m,\sigma^2/2\mu)$.

The operator runs either mode 0 or mode 1. Running mode 1 instead of mode 0 earns the flow $\Phi(Q)$ per unit time, for example $\Phi(Q)=C_0(Q)-C_1(Q)$ when the modes have operating costs $C_0,C_1$. Switching 0 → 1 costs $F^+\ge 0$, switching 1 → 0 costs $F^-\ge0$, and future payoffs are discounted at rate $\rho>0$. Assume $\Phi$ is increasing with a single root $\hat Q$, so mode 1 is the high-demand mode. If $\Phi$ decreases instead, the app (`index.html`) mirrors the problem through $Q\mapsto 2m-Q$; this is valid because the OU process is symmetric about $m$.

## 2. Bellman equation and the continuation regions

Let $V_0(Q)$ and $V_1(Q)$ be the optimal values when currently in mode 0 or mode 1, measured relative to staying in mode 0 forever. The optimal policy has two triggers $Q_L<Q_H$: switch 0 → 1 when $Q$ reaches $Q_H$, and 1 → 0 when it reaches $Q_L$. Between switches, Itô's lemma applied to the Bellman equation gives

$$
\tfrac12\sigma^2V_0''+\mu(m-Q)V_0'-\rho V_0=0,\qquad Q<Q_H,
$$

$$
\tfrac12\sigma^2V_1''+\mu(m-Q)V_1'-\rho V_1+\Phi(Q)=0,\qquad Q>Q_L.
$$

The diffusion coefficient is $\tfrac12\sigma^2$, not $\tfrac12\sigma^2Q^2$: this is the only change to the operator compared with the proportional-volatility process. It changes which special functions solve the equation.

## 3. Homogeneous solutions: parabolic cylinder functions

With $z=\sqrt{2\mu}\,(Q-m)/\sigma$ and $\nu=\rho/\mu$, the homogeneous equation becomes $V_{zz}-zV_z-\nu V=0$. Writing $V=e^{z^2/4}u(z)$ gives Weber's equation $u''+\big(\tfrac12-\nu-\tfrac{z^2}{4}\big)u=0$, whose solutions are the parabolic cylinder functions $D_{-\nu}(\pm z)$. Two convenient positive solutions are

$$
G_+(Q)=\Gamma(\nu)\,e^{z^2/4}D_{-\nu}(-z)=\int_0^\infty t^{\nu-1}e^{zt-t^2/2}\,dt,\qquad G_-(Q)=G_+\big|_{z\to -z}.
$$

$G_+$ is increasing and tends to 0 as $Q\to-\infty$; $G_-$ is decreasing and tends to 0 as $Q\to+\infty$. Their derivatives are $G_\pm'(Q)=\pm\frac{\sqrt{2\mu}}{\sigma}\int_0^\infty t^{\nu}e^{\pm zt-t^2/2}dt$. Equivalently, they are combinations of $M\!\big(\tfrac{\nu}{2},\tfrac12,\tfrac{z^2}{2}\big)$ and $z\,M\!\big(\tfrac{\nu+1}{2},\tfrac32,\tfrac{z^2}{2}\big)$, Kummer functions in $z^2$ rather than in $1/Q$.

## 4. Boundary conditions and the general solution

- As $Q\to-\infty$, the option to switch into mode 1 becomes worthless, so $V_0$ cannot contain the exploding $G_-$: $V_0(Q)=A\,G_+(Q)$.

- As $Q\to+\infty$, the option to leave mode 1 becomes worthless, so $V_1$ approaches the value of staying forever: $V_1(Q)=B\,G_-(Q)+P(Q)$.

The particular solution is the expected discounted payoff of never leaving mode 1. By the Feynman–Kac formula,

$$
P(Q)=\mathbb E\!\left[\int_0^\infty e^{-\rho s}\,\Phi(Q_s)\,ds\ \Big|\ Q_0=Q\right]=\int_0^\infty e^{-\rho s}\,\mathbb E\big[\Phi(Q_s)\big]\,ds,
$$

using the Gaussian transition law above. This is the expectation of $\Phi$, not $\Phi$ of the expectation. For polynomial payoffs it has a closed form. For $\Phi=a+bQ+cQ^2$,

$$
P(Q)=\frac{a+bm+c\,(m^2+\tfrac{\sigma^2}{2\mu})}{\rho}+\frac{(b+2cm)(Q-m)}{\rho+\mu}+\frac{c\,\big((Q-m)^2-\tfrac{\sigma^2}{2\mu}\big)}{\rho+2\mu}.
$$

For general $\Phi$, the app evaluates the double integral numerically, with Gauss–Hermite quadrature over the Gaussian (details in Section 5).

## 5. Value matching and smooth pasting

At the triggers, the value of switching (net of its cost) must equal the value of staying, and the slopes must match:

$$
A\,G_+(Q_H)=B\,G_-(Q_H)+P(Q_H)-F^+,\qquad B\,G_-(Q_L)+P(Q_L)=A\,G_+(Q_L)-F^-,
$$

$$
A\,G_+'(Q_H)=B\,G_-'(Q_H)+P'(Q_H),\qquad A\,G_+'(Q_L)=B\,G_-'(Q_L)+P'(Q_L).
$$

The four equations in $(A,B,Q_L,Q_H)$ are linear in $A,B$; eliminating them leaves two nonlinear equations in $(Q_L,Q_H)$.

**How the app solves them.** Inside the band $(Q_L,Q_H)$ both value functions satisfy their equations, so the difference $D=V_1-V_0$ satisfies

$$\tfrac12\sigma^2D''+\mu(m-Q)D'-\rho D+\Phi(Q)=0,\qquad Q_L<Q<Q_H,$$

with $D(Q_H)=F^+$, $D(Q_L)=-F^-$ (value matching) and $D'(Q_H)=D'(Q_L)=0$ (smooth pasting, since $V_0'=V_1'$ at both triggers). For trial triggers, the app integrates this ODE inward from each end with fourth-order Runge–Kutta and requires the two solutions to meet with equal value and slope at $\hat Q$. Newton's method on that mismatch, with a continuation that scales the costs up from near zero, gives $(Q_L,Q_H)$. The coefficients $A,B$ then follow from smooth pasting. Shooting inside the band avoids evaluating $G_\pm$ far into their exponential tails, which keeps the solve stable for small discount rates and wide bands. When one trigger lies beyond ±8 long-run SD, the app solves the remaining one-sided problem with $G_\pm$ and $P$ directly.

$P(Q)$ is computed as $\bar E/\rho+R(Q)$, where $\bar E$ is the long-run mean of $\Phi$ and $R$ is the integral of the transient part $e^{-\rho s}\big(\mathbb E[\Phi(Q_s)]-\bar E\big)$ over graded Gauss–Legendre panels in $s$. This stays accurate when $\rho\ll\mu$.

## 6. Properties used as checks

- **No switching costs.** As $F^\pm\to0$, both triggers converge to the break-even point $\hat Q$ with $\Phi(\hat Q)=0$. Subtracting the two equations at a common trigger with $V_0=V_1$, $V_0'=V_1'$ and $V_0''=V_1''$ gives $\Phi(\hat Q)=0$. With free switching, the myopic rule is optimal.

- **Switching regions.** Optimality requires $\Phi(Q)\ge\rho F^+$ wherever the policy enters mode 1 and $\Phi(Q)\le-\rho F^-$ wherever it leaves, so $Q_L<\hat Q<Q_H$. The app also verifies on a grid that neither value function can be improved by switching early.

- **One-sided solutions.** If one switch is never worth making within ±8 long-run standard deviations of $m$, the corresponding coefficient is zero and the app solves for the single remaining trigger.

## 7. Estimating the OU process from data

The exact discretization over a step $\Delta_i$ is

$$
Q_{i}=m+(Q_{i-1}-m)e^{-\mu\Delta_i}+\varepsilon_i,\qquad \varepsilon_i\sim N\!\Big(0,\tfrac{\sigma^2}{2\mu}\big(1-e^{-2\mu\Delta_i}\big)\Big).
$$

With equal spacing $\Delta$, this is an AR(1) model $Q_i=a+bQ_{i-1}+\varepsilon_i$, and the conditional maximum-likelihood estimates are $\hat\mu=-\ln\hat b/\Delta$, $\hat m=\hat a/(1-\hat b)$ and $\hat\sigma^2=2\hat\mu\,\hat s^2/(1-\hat b^2)$, where $\hat s^2$ is the residual variance. With irregular spacing, the app maximizes the exact Gaussian likelihood numerically. Standard errors come from the observed information (a numerical Hessian). Mean reversion requires $0<\hat b<1$. The speed $\mu$ is weakly identified unless the sample spans many half-lives $\ln2/\mu$.

## 8. Differences from Guo, Chow & Schonfeld (2018)

- The paper's valuation equations use $\tfrac12\sigma^2Q^2V''$ (the process $dQ=\mu(m-Q)dt+\sigma Q\,dW$), which leads to Kummer functions in $1/Q$. With $\sigma=7$, that process has infinite variance and Eq. (41) explodes when simulated. The simulated paths in the paper behave like additive OU, which is the model solved here.

- Here $P(Q)$ is the exact expectation $\mathbb E[\Phi(Q_s)]$ under the Gaussian transition law, rather than $\Phi(\mathbb E[Q_s])$.

- The boundary conditions are imposed at $Q\to\pm\infty$, where each option becomes worthless. With zero switching costs, the triggers coincide with the root of $\Phi$.

- The smooth-pasting conditions use $\tfrac{d}{dQ}$ of the full solution. In the paper's Eq. (38), the derivative rows are missing a factor $1/Q$ and have the wrong sign on $\partial P/\partial Q$.

## References

- Dixit, A. (1989). Entry and exit decisions under uncertainty. *Journal of Political Economy* 97(3), 620–638.

- Dixit, A., Pindyck, R. (1994). *Investment under Uncertainty*. Princeton University Press.

- Guo, Q.-W., Chow, J.Y.J., Schonfeld, P. (2018). Stochastic dynamic switching in fixed and flexible transit services as market entry-exit real options. *Transportation Research Part C: Emerging Technologies* 94, 288–306. https://doi.org/10.1016/j.trc.2017.08.008

- Abramowitz, M., Stegun, I. (1964). *Handbook of Mathematical Functions*, ch. 19 (parabolic cylinder functions).
