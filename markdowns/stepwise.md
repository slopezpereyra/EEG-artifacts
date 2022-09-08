
## Stepwise analysis algorithm
### Justification
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
