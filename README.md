# Adaptive Truncation Table Generator

- `com_poisson_iterations.R`, has dependecy:
  - `aux/aux.R`: function `log_diff_exp`;
  - `stan/comp_Z_bp.stan`: Bounding pairs approach for the COM-Poisson normalising constant;
  - `stan/comp_Z_brute.stan`: The COM-Poisson normalising constant by brute force;
  - `stan/comp_Z_brms_fixed.stan`: The COM-Poisson normalising constant by brms with index error corrected and M increased from 10^4 to 10^6;

![image](https://github.com/user-attachments/assets/f322c6ab-8290-478e-be21-1095ef6d2854)

- `com_poisson_error.R`, has dependency:
  - `aux/aux.R`: function `log_diff_exp`;
  - `stan/comp_Z_bp.stan`: Bounding pairs approach for the COM-Poisson normalising constant;
  - `stan/comp_Z_brute.stan`: The COM-Poisson normalising constant by brute force;

![image](https://github.com/user-attachments/assets/102da055-8e66-4474-bb13-7571f686965e)


- `domination_bp_over_stt.R`, has dependency:
  - `domination_bp_over_stt.stan`: Implementation of Bounding pairs and Sum-to-threshold for the problem;

![image](https://github.com/user-attachments/assets/04113b96-cc48-4d0b-b0f6-c29b9a3098fd)
