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

![image](https://github.com/user-attachments/assets/d6cdf301-8cd9-46f4-b7cf-b25ab33a637c)

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
![image](https://github.com/user-attachments/assets/abc0e82d-6409-42d2-ac67-6e5faf0aae63)

- `domination_bp_over_stt.py`, has dependency:
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sum_to_threshold_mp.py`: The Sum-to-threshold approach with mpmath structure;

Run:

```bash
python domination_bp_over_stt.py
```

![image](https://github.com/user-attachments/assets/990b1f49-02a7-4091-a415-e48686614627)

