from tabulate import tabulate

# Data for Group A and Group B
group_a = [
    ['A1', 'A2'],
    ['A3', 'A4'],
    ['A5', 'A6']
]

group_b = [
    ['B1', 'B2'],
    ['B3', 'B4'],
    ['B5', 'B6']
]

# Combine group A and group B data with index
data = []
for idx, (a, b) in enumerate(zip(group_a, group_b), start=1):
    data.append([idx, a[0], a[1], b[0], b[1]])

# Headers for the table
headers = ["Index", "A Col 1", "A Col 2", "B Col 1", "B Col 2"]

# Print table with tabulate in a 'grid' format
print(tabulate(data, headers, tablefmt="fancy_grid"))