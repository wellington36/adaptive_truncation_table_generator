# Adaptive Truncation Table Generator

- `comp_Z_iterations.py`, has dependency:
  - `utils/utils.py`: function `logdiffexp`;
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;
  - `src/brute_mp.py`: The Brute force approach with mpmath structure
  - `stan/comp_Z_brms_fixed.stan`: The COM-Poisson normalising constant by brms with index error corrected and M increased from _10^4_ to _10^6_;

Run (check if auxiliary libraries are installed):
```bash
python comp_Z_iterations.py
```

![image](https://github.com/user-attachments/assets/0455d6f9-36d7-4713-80b4-26bfb2e95ef9)


- `comp_Z_errors.py`, has dependency:
  - `utils/utils.py`: function `logdiffexp`;
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;
  - `src/brute_mp.py`: The Brute force approach with mpmath structure
  - `stan/comp_Z_brms_fixed.stan`: The COM-Poisson normalising constant by brms with index error corrected and M increased from _10^4_ to _10^6_;

Run (check if auxiliary libraries are installed):
```bash
python comp_Z_errors.py
```
![image](https://github.com/user-attachments/assets/1c944657-f6f2-4221-a7ef-16af12c7c9db)



- `domination_bp_over_stt.py`, has dependency:
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sum_to_threshold_mp.py`: The Sum-to-threshold approach with mpmath structure;

Run:

```bash
python domination_bp_over_stt.py
```

![image](https://github.com/user-attachments/assets/b3a2b72a-5471-40ce-8ab2-82f271a20f0b)

- `p_series_iterations.py`, has dependency:
  - `src/integration_bound_mp.py`: The Bounding pairs with integral test approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;

Run:

```bash
python p_series_iterations.py
```

![image](https://github.com/user-attachments/assets/50331758-2663-411b-83ca-ed1f5fcdd2ef)

