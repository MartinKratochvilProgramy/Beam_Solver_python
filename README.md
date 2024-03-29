# How to use
Install requirements by running:
```
pip install -r requirements.txt
```
Run the program by:
```
python main.py
```

# User instructions:

Run by opening and running window.py - after running, window with user input appears. In the user input window, add lenght L, support coordinates Xa and Xb, moment of inertia Jz and apply forces. Note that force or linear load facing downwards should be negative. After hitting the Solve! button, program computes numerically shear forces and bending moments in the beam. Note that number of divisions is set at 10 000. Also since the calculation is numerical and to solve displacements initial deflection and finite deflection step are defined, user has to input somewhat reasonable values, otherwise the calculation might either take too long or likely not converge at all.

![image](https://user-images.githubusercontent.com/94861828/148649023-289b3766-a5f5-4034-b3ea-6193af2bb383.png)

# Results:
Result is a matplotlib graph with input forces and another three graphs with shear forces, bending moments and displacements

![image](https://user-images.githubusercontent.com/94861828/148651299-905c2f02-bc02-4626-9fcb-ba231a969bb1.png)

# Disclaimer!
Program is not inteded for commercial use, I will not be responsible for any potential failures caused by the use.





