# Derivations
This file contains some theoretical derivations of the collision equations used in this simulation, in case you would like to see those.
That is, we won't look at general mathematical derivations of rotation matrices or what not, but at the (much more interesting) physical derivations.

#### 1. Momentum transfer in perfectly elastic collisions ($e=1$) in 1D.

As in the simulation, we consider one object to be moving while the other is stationary. Here we will not also consider the effect of the gravitational force though. Say object 1 has velocity $v_1$ (and $v_2=0$) before the collision, and the objects have masses $m_1$ and $m_2$. We want to calculate the velocities of both objects after the collision, $v'_1$ and $v'_2$. 

This can be done by considering momentum and energy conservation. Momentum conservation gives

$$ m_1 v_1 = m_1 v'_1 + m_2 v'_2, $$

while energy conservation gives

$$ \frac{1}{2} m_1 v_1^2 = \frac{1}{2} m_1 (v'_1)^2 + \frac{1}{2} m_2 (v'_2)^2. $$

To combine these, we rewrite the momentum equation as 

$$ v'_2 = \frac{m_1}{m_2} (v_1 - v'_1) $$

and substitute this in the energy equation to find

$$ \frac{1}{2} m_1 \left(v_1^2 - (v'_1)^2 \right) = \frac{1}{2} m_2 (v'_2)^2 = \frac{m_1^2}{2m_2} (v_1 - v'_1)^2 $$

so that

$$ v_1^2 - (v'_1)^2 = \frac{m_1}{m_2} ( v_1^2 + (v'_1)^2 - 2v_1v'_1 ) $$

which can be rewritten as

$$ \left( \frac{m_1}{m_2} + 1 \right) (v'_1)^2 - 2 \frac{m_1}{m_2} v_1 v'_1 + \left( \frac{m_1}{m_2} - 1 \right) v_1^2 = 0. $$

This can be solved for $v'_1$ using the ABC-formula, giving

$$ v'_1 = \frac{2\frac{m_1}{m_2}v_1 \pm \sqrt{ 4\frac{m_1^2}{m_2^2}v_1^2 - 4\left( \frac{m_1}{m_2} + 1 \right)\left( \frac{m_1}{m_2} - 1 \right)v_1^2 }}{2\left( \frac{m_1}{m_2} + 1 \right)} \\ \: \\
= \frac{\frac{m_1}{m_2} \pm \sqrt{\frac{m_1^2}{m_2^2} - \left( \frac{m_1^2}{m_2^2} - 1\right)}}{ \frac{m_1}{m_2} + 1 } v_1 \\ \: \\
= \frac{m_1 \pm m_2}{m_1+m_2} v_1 .$$

From this follows $v'_2$ as well:

$$ v'_2 = \frac{m_1}{m_2} (v_1 - v'_1) = \frac{m_1}{m_2} (v_1 - \frac{m_1 \pm m_2}{m_1+m_2} v_1) \\ \: \\ = \frac{m_1}{m_2} \left( \frac{m_2 \mp m_2}{m_1+m_2} \right) v_1 \\ \: \\ = \frac{m_1 \mp m_1}{m_1 + m_2} v_1 .$$

We have found two solutions here. Note that for the 'upper' one (using the $+$ in $v'_1$ and $-$ in $v'_2$), we get $v'_1 = v_1$ and $v'_2 = 0$: this is what happens when a collision does not actually occur, meaning the velocities before and after remain the same. This is a valid solution, but not the one we are interested in. So after a collision, the velocities are given by the other 'lower' solution:

$$ v'_1 = \frac{m_1 - m_2}{m_1+m_2} v_1, \:\:\:\:\: v'_2 = \frac{2m_1}{m_1 + m_2} v_1 .$$

Note how in the case of equal masses $m_1=m_2$, we find $v'_1=0$ and $v'_2=v_1$, as one would expect. 

$$ \: $$

#### 2. Momentum transfer in inelastic collisions in 1D.

We now consider inelastic collisions, meaning the energy conservation from before no longer holds, as the kinetic energy leaks away into other forms of energy. We do still have momentum conservation though. To understand why, let us e.g. consider conversion of kinetic energy into thermal energy. This thermal energy would increase the velocities and thus momenta of the particles making up the object. However, given the nature of thermal movement, the net momentum of all particles remains the same (zero, or whatever momentum the object as a whole is moving with). So while the momenta within all increase, they cancel each other out in every direction; this is possible due to the vector nature of momentum. So the total momentum of internal particles stays the same, meaning the momenta of the colliding objects will remain constant; the momentum cannot leak away through thermal energy changes. Generally speaking, the fact that momentum is a vector also prevents it from leaking away when energy is converted into something else. (Also, for the case of single, elemental particles colliding, there is of course e.g. no temperature associated with them, and so energy cannot leak away there either; particle collisions always have to be perfectly elastic.)

Instead of energy conservation, we add the definition of the coefficient of restitution as one of our equations:

$$ e = \frac{|v'_1 - v'_2|}{|v_1 - v_2|} .$$

Since we consider the case of $v_2=0$, this can be written as

$$ v'_1 - v'_2 = \pm e v_1 ,$$

where the $\pm$ has to be chosen such that both sides have the same sign.
Now combining this with momentum conservation, we find

$$ v'_1 = \pm e v_1 + v'_2 = \pm e v_1 + \frac{m_1}{m_2}(v_1 - v'_1) $$

so 

$$ v'_1 \left( \frac{m_1}{m_2} + 1\right) = v_1 \left( \frac{m_1}{m_2} \pm e \right) $$

$$ \Rightarrow v'_1 = \frac{\frac{m_1}{m_2} \pm e}{\frac{m_1}{m_2} + 1} v_1 = \frac{m_1 \pm e m_2}{m_1 + m_2} v_1 .$$

And then 

$$ v'_2 = v'_1 \mp e v_1 = \frac{m_1 \pm e m_2}{m_1 + m_2} v_1 \mp \frac{e(m_1+m_2)}{m_1+m_2} v_1 \\ \: \\
= \frac{m_1 \mp e m_1}{m_1+m_2} v_1 .$$

By comparing this to the perfectly elastic case $e=1$, the 'upper' solution is again valid if no collision actually occurs, so we want the 'lower' solution to find the velocities after a collision:

$$ v'_1 = \frac{m_1 - e m_2}{m_1+m_2} v_1, \:\:\:\:\: v'_2 = \frac{(1+e)m_1}{m_1 + m_2} v_1 .$$