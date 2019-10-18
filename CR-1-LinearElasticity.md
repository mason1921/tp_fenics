# Report 1 - MU5MES01 - 2019/20 
*S.Brisard, D.Duhamel, C.Maurini*

The key goals of the first part of the class is to 
 - Be able to implement a finite element solver for linear elasticity.
 - Be able to perform a convergence analysis.
    
Your report should summarize and present synthetically your work on these items. We give below some hints on how to write the report. Personalized analyses and comments are particularly welcome. Your are not obliged to follow the following format step and by step. But you should include in your report the key concepts and results. 



**Some suggestions:**
 - Write concisely and effectively.
 - Comment your results.
 - The quality of the figures is important.
 - Report only the minimal number of figures (of excellent quality) to effectively communicate your results.
 - You can write in English or French.
 - Use Latex for writing your report. 


**Important informations:**
  - Deadline: For the final version **Thursday 24 october, 23h59**. We strongly suggest to have a preliminary version of the codes availble for Friday 18 october. 
  - **The maximal length of the report is 4 pages. **
  - To submit your report: 
      - a paper copy of your report should be returned by Thursday 24/10 at 8:30. One paper version per group is sufficient.
      - An electronic version should be submitted to github. Proceed as follows to create the work and submission repository for your group:
        - **Only one of the two students** of your group will go to  https://classroom.github.com/g/XBKAl9Tu, accept the assignement and create a `new team`, naming the team as `NAMESTUDENT1-NAMESTUDENT2`.
        - **Once the first student has create the team**, the second student goes to https://classroom.github.com/g/XBKAl9Tu, accepts the assignement and asks to join the team with his name (do not create another team, there should one team for group).
      - In your "group" repository you should: 
          1. Create a directory called `CR1`
          2. Put your report in the pdf form named as `MES01-CR1-studentname1-studentname2.pdf` (file with a different naming scheme will not be accepted and evaluated. 
          3. Put all your files you used to obtain your results in `CR1/src` (namely the *.py and *.ipynb files)
  - We will evaluate the quality of the presentation (language, typesetting, and figures). Being able to effectively communicate your results is important for your future.
  - We require you to be able to use git at least to push your data to the repository. This is the main reason why we ask to submit your report on the github platform. We will not accept submissions by mail.
      
      


# Linear elasticity and convergence analysis

## Part I - Cantilever beam

We consider the problem presented in the notebook `1.0-Cantilever_beam.ipynb`

Your report should include:
- A variational formulation of the problem.
- The statement of the type of discretisation used to solve the problem (type of finite elements).
- A figure with a comparison between the analytical solution for the shear force on the beam when using Lagrange elements of degree 1 and 2, with a tentative discussion of the results.
- A figure reporting the convergence analysis for linear and quadratic elements. You can follow as guidelines the fully commented converged analysis of Section 5.5 of the FEniCS tutorial, see in particular Section 5.5.4. Comment the results on the convergence rate observed in your numerical experiments. You can use the class `Cantilever` we provided to perform the analysis. 
- Comments and conclusions (the `Questions` in the notebook are possible hints to produce pertinent comments and remarks)

Any pertinent complement to the previous analysis will be appreciated in the grading and will be the object of discussion at the final oral examination. 

**Ideas for advanced study for part I**

Perform this part only if you feel your level sufficient advanced; this part is not required if you are aiming at grades < 14/20.

- Perform the convergence analysis for different values of the Poisson ratio. Report the convergence plot for $\nu=0.499999999$ and comment the results.


## Part II -  Cook membrane
Write a symple `cook.py` script solving the 
the Cook's membrane. In this classical test case, we set $E = 250, ν = 0.4999, V = 100.$
You can use the script `linear_elasticity.py` we provided as a starting point.
To generate the mesh you can use the example below (modify the points). 
Produce a figure illustrating the results for the displacement and the von Mises stress.
```python
domain_vertices = [dolfin.Point(0.0, 0.0125),
                  dolfin.Point(0.0, -0.025),
                  dolfin.Point(1.0, -1.0),
                  dolfin.Point(2.0, 0.0),
                  dolfin.Point(1.0, 1.0)]
domain = mshr.Polygon(domain_vertices)
ndiv = 30
mesh = mshr.generate_mesh(domain,ndiv)
```

<img src="images/cook.png" alt="Drawing" style="width: 300px;"/>


# References
[1] H.P.Langtangen, A.Logg, Solving PDEs in Minutes - The FEniCS Tutorial Volume I, Springer 2016 https://www.springer.com/gp/book/9783319524610

[2] B.Szabó, I.Babuška, Introduction to Finite Element Analysis: Formulation, Verification and Validation, Wiley 2011, ISBN: 978-0-470-97728-6


```python

```
