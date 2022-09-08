
## Stepwise analysis algorithm

### Why stepwise analysis?

A valid question is that which interrogates about the need for stepwise analysis when one could, in principle, perform a single analysis on the whole EEG record. The answer to this question lies in how CAPA, the statistical procedure implemented for anomaly detection, works.

The starting point of CAPA analysis is to produce an estimation $\hat{\theta}$ of the parameters $\theta$ of the standard (or not anomalous) data. Such estimation is then used to infer epidemic changes in the distribution via minimization of the cost function 

$$ \sum \limits_{t\notin \cup \left[{\tilde{s}}_i+1,{\tilde{e}}_i\right]}\mathcal{C}\left({\mathbf{x}}_t,{\hat{\theta}}_0\right)+\sum \limits_{j=1}^{\hat{K}}\left[{\min}_{{\tilde{\theta}}_j}\left(\sum \limits_{t={\tilde{s}}_j+1}^{{\tilde{e}}_j}\mathcal{C}\left({\mathbf{x}}_t,{\tilde{\theta}}_j\right)\right)+\beta \right], $$

where $e_i - s_i \geq \hat{l}$, with $\hat{l}$ the minimum anomalous segment length, $\mathcal{C}$ some cost function and $\beta$ an adecuate penalty. (For more information, see the original paper by [Fisch, Eckley and Fearnhead](https://onlinelibrary.wiley.com/doi/full/10.1002/sam.11586).)

Because this cost function depends on the estimated parameters $\hat{\theta}$, and because sleep data is subject to regular shifts in the real distribution parameters (e.g., in the transition from one sleep stage to the other), it is best to estimate the parameters of more or less brief regions of the distribution locally. In fact, this is how artifact detection is carried out by human experts: by detecting deviations from relatively local patterns (e.g., a thirty seconds epoch).

### Algorithm logic
Let $t$ be the duration in seconds of a full record, and $s$ the duration in seconds of a single step. Let also $T(a)$ return the duration in seconds of a subset $a$ of the record. Then

$i)$ If $s | t$, $q$ analysis will cover the full record with $t = s \times q$ and a sequence of analyzed subsets $a_0, a_1, ...., a_q$, with $T(a_i) = s$  for all $i$.

$ii)$ If $s \nmid t$, $t = s \times q + r, 0 \leq r < q$ with a sequence $a_1, a_2, ..., a_q, a_{q+1}$ of steps, with $T(a_i) = s$ for all $i \neq q + 1$, and $T(a_{q + 1}) = r$.

From this reasoning sprouts the stepwise algorithm. It is designed to perform $q$ analysis of length $s$ if $s | t$, or $q$ analysis of length $s$ and an extra one of length $r$ if $s \nmid t$. Thus, it analyizes the whole EEG record. 

### Pseudo-code implementation

```python
--- ARGS

eeg = an eeg dataset
s = duration of each step of analysis

--- FUNC
def stepwise_analysis(eeg, s):
	# Start end seconds of a_0 (first sequence)
	s = eeg.time[1] # First second of the provided EEG
	e = start + s
	
	last_second = eeg.time[-1] # Last second of the provided EEG
	t = last_second - s # Duration of the EEG
	steps = t %/% s # quotient
	r = t %% s # remainder
	
	if r != 0:
		steps += 1

	for i in (1, 2, ..., steps):
		if (e > last_second): # Logically, only occurs on the last step if 
							  # r != 0 on the last step
			e = s + r 		  # Make (q + 1)th analysis of length r
			 			      # (not s).
			
		analyze(eeg, s, e) # Analyze the EEG from second s to e
		s = e
		e += s
```
