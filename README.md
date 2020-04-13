# MATF-RI20-PIGA-za-Dizajn-Radio-Mreze

# :tophat: Project Title

Ovaj rad je baziran na realnom problemu kombinatorne optimizacije, kao primer koji pokazuje kako genetski algoritam (GA) može biti paralelizovan na efikasan način. Problem koji ćemo razmatrati ima za cilj odredivanje najboljeg skupa lokacija radio predajnika kako bi se pokrio što veći prostor po optimalnoj ceni. Dati problem ćemo prvo testirati na Standardnom Genetskom Algoritmu (SGA), a potom na Ostrvskom modelu Genetskog Algoritma (IGA) koji ćemo paralelizovati na višejezgarnom procesoru.

## :floppy_disk: Prerequisites:
```
OpenMP library
```
```
OpenGL library
```
```
C++
```

## :wrench: Setup: 

`$ OpenMP library - sudo apt-get install libomp-dev`

`$ OpenGL library - sudo apt-get install mesa-utils`


## :heavy_check_mark: Startup:
`$ g++ std=c++14 sga.cpp -o sga -lGL -lGLU -lglut -fopenmp`

`$ g++ std=c++14 island.cpp -o island -lGL -lGLU -lglut -fopenmp`


## :mortar_board: Autors:
* Tijana Todorov       | e-mail &bull; tijana.todorov710@gmail.com
* David Nedeljkovic    | e-mail &bull; dnedeljkovic710@gmail.com
