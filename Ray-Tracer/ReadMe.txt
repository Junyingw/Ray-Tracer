Assignment #3: Ray tracing

==== BUILD ON MAC OS X ====

> unzip 420_assignment3.zip
> cd pic 
> make 
> cd ..
> cd assign3
> make
> ./assign3 <screenfile>[1 for recursive reflection][1 for antialiasing][jpegname]

============== PAY ATTENTION ============

1) You can set the times of recursive reflection via source code manually, the default recursive reflection level is 3 
 
2) Several ways of running the program:

>./assign3 <screenfile> 
(Just the still image, no antialiasing, no recursive reflection )

>./assign3 <screenfile>[1 for recursive reflection][1 for antialiasing]
(You can set antialiasing and recursive reflection. i.e ./assign3 screenfile.txt 1 1 )

>./assign3 <screenfile>[1 for recursive reflection][1 for antialiasing][jpegname]
(You can set antialiasing, recursive reflection as well as saving the jpeg file  i.e ./assign3 test1.txt 1 1 0001.jpeg )

3) Keyboard: 'q' and 'Q' for easy quit the program
 
=================== MANDATORY FEATURES ========================

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                  yes

2) Ray tracing sphere                     yes

3) Triangle Phong Shading                 yes

4) Sphere Phong Shading                   yes

5) Shadows rays                           yes

6) Still images                           yes
   
7) Extra Credit (up to 20 points)

1. Recursive reflection (Default value is 3, you can set via source code)

2. Good antialiasing (Use average filter to reduce jags, may cause image blur)

===================== SCENE FILE OUTPUT ==============================

Name          Test Scene      Recursive Reflection   Antialiasing
-----         -----------     --------------------   -------------
0001	      test1.scene		no                   no
0002	      test1.scene		yes		     no
0003	      test1.scene		no                   yes
0004	      test1.scene		yes                  yes

0005	      test2.scene		no                   no
0006	      test2.scene		yes		     no
0007	      test2.scene		no                   yes
0008	      test2.scene		yes                  yes             

0009	      spheres.scene		no                   no
0010	      spheres.scene		yes		     no
0011	      spheres.scene		no                   yes
0012	      spheres.scene		yes                  yes

0013	      table.scene		no                   no
0014	      table.scene		yes		     no
0015	      table.scene		no                   yes
0016	      table.scene		yes                  yes             

0017	      siggraph.scene		no                   no
0018	      siggraph.scene		yes		     no
0019	      siggraph.scene		no                   yes
0020	      siggraph.scene		yes                  yes             




