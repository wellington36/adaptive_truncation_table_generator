# Adaptive Truncation Table Generator

- `comp_Z_iterations.py`, has dependency:
  - `utils/utils.py`: function `logdiffexp`;
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;
  - `src/fixed_mp.py`: Evaluate a fixed numbers of terms
  - `stan/comp_Z_brms_fixed.stan`: The COM-Poisson normalising constant by brms with index error corrected and M increased from _10^4_ to _10^6_;

Run (check if auxiliary libraries are installed):
```bash
python comp_Z_iterations.py
```

![image](https://github.com/user-attachments/assets/4e9c581d-a9f5-4d73-95f7-874b2daa066f)


- `comp_Z_errors.py`, has dependency:
  - `utils/utils.py`: function `logdiffexp`;
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;
  - `src/fixed_mp.py`: Evaluate a fixed numbers of terms
  - `stan/comp_Z_brms_fixed.stan`: The COM-Poisson normalising constant by brms with index error corrected and M increased from _10^4_ to _10^6_;

Run (check if auxiliary libraries are installed):
```bash
python comp_Z_errors.py
```
![image](https://github.com/user-attachments/assets/74f3cb00-c63f-4f90-b81f-ca4fcc0976d9)



- `domination_bp_over_stt.py`, has dependency:
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/sum_to_threshold_mp.py`: The Sum-to-threshold approach with mpmath structure;

Run:

```bash
python domination_bp_over_stt.py
```

![image](https://github.com/user-attachments/assets/8477a9e9-c1b0-4390-b374-350858f3fc6c)

- `comp_Z_approx_vs_bounding.py`, has dependency:
  - `src/bounding_pairs_mp.py`: The Bounding pairs approach with mpmath structure;
  - `src/fixed_mp.py`: Evaluate a fixed numbers of terms

Run:

```bash
python comp_Z_approx_vs_bounding.py
```

![image](https://github.com/user-attachments/assets/a1814353-ccf8-4ece-8a24-925b344af599)

- `p_series_iterations.py`, has dependency:
  - `src/integration_bound_mp.py`: The Bounding pairs with integral test approach with mpmath structure;
  - `src/sequential_mp.py`: The Sequential approach with brute force and mpmath structure;


Run:

```bash
python p_series_iterations.py
```

![image](https://github.com/user-attachments/assets/50331758-2663-411b-83ca-ed1f5fcdd2ef)

- Table of mcmc (made in https://github.com/wellington36/MCMC_COMPoisson)

![image](https://github.com/user-attachments/assets/d973cc3f-673c-409f-90ef-0765b5dc5191)
