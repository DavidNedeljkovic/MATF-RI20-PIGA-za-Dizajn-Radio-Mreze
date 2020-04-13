#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <map>
#include <tuple>
#include <set>
#include <cmath>
#include <omp.h>
#include <GL/glut.h>

const int MAX_ITERATIONS = 251;
const int POPULATION_SIZE = 160;
const int ELITE_SIZE = 24;
const int MUTATION_PROB = 7;
const int INIT_PROBABILITY = 35;
const int NEGATIVE_INFINITY = -1E9;
const int TOURNAMENT_SIZE = 6;
const double OPSEG = 20.5;

static void on_keyboard(unsigned char key, int x, int y);
static void on_reshape(int width, int height);
static void on_display(void);

static void init(void)
{
    // Obavlja se OpenGL inicijalizacija.
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0.0, 0.0, 0.0, 0.0);
}


class Jedinka
{
public:

	std::vector<bool> kod_jedinke;
	std::vector<std::tuple<double,double>> korisceni_releji;
	std::map<std::tuple<int,int>,int> mapa_lokacija;
	double fitnes;

	Jedinka(std::map<std::tuple<double,double>,std::vector<std::tuple<int,int>>> mapa_releja,
		 std::vector<std::tuple<double,double>> releji,
		 int pokrivenost)
	{
		// kreiranje koda jedinki
		kod_jedinke.resize(releji.size());
		for (int i = 0; i < releji.size(); i++) {
			bool tmp = ((std::rand() % 100) > INIT_PROBABILITY);
		    kod_jedinke.at(i) = tmp;
		}
		popraviKodJedinke(releji);
		fitnes = fitnesFunkcija(pokrivenost, mapa_releja, releji);
	}
    // funkcija koja proverava da li postoji bar jedna 1
	// ako ne onda ubacuje jednu 1 na random poziciju
	void popraviKodJedinke(std::vector<std::tuple<double,double>> releji)
	{
		for (int i = 0; i < releji.size(); i++) {
		    if (kod_jedinke.at(i))
		    	return;
		}
	    int i = std::rand() % releji.size();
		kod_jedinke.at(i) = true;
	}
	// funkcija koja racuna fitnes vrednost za svaku jedinku
	double fitnesFunkcija(int pokrivenost,
		std::map<std::tuple<double,double>,std::vector<std::tuple<int,int>>> mapa_releja,
		std::vector<std::tuple<double,double>> releji)
	{
		mapa_lokacija.clear();
		korisceni_releji.clear();

		//kreiranje liste koriscenih releja
		for(int i=0; i < kod_jedinke.size(); i++){
			if(kod_jedinke.at(i)){
		    	korisceni_releji.push_back(releji.at(i));
		    }
		}
		//kreiranje mape lokacija
		for(int i=0; i < korisceni_releji.size(); i++){
		   std::vector<std::tuple<int,int>> ml = mapa_releja.at(korisceni_releji.at(i));
		   for(int j=0; j < ml.size(); j++){
		     	if(mapa_lokacija.find(ml.at(j)) == mapa_lokacija.end()){
					mapa_lokacija.insert({ ml.at(j), 1 });
		        }
				else{
					mapa_lokacija.at(ml.at(j)) += 1;
			 	}
			}
		}
		//racunanje fitnesa
		double kvadrat  = pow(100.0*(mapa_lokacija.size())/(pokrivenost*1.0), 2);
		double rezultat = (kvadrat / (korisceni_releji.size()*1.0));
		return rezultat;
	}
};

std::vector<Jedinka> populacija;
std::vector<Jedinka> nova_populacija;

// funkcija koja vrsi turnirsku selekciju i vraca pobednika
int selekcija()
{
	double max = NEGATIVE_INFINITY;
	int pobednik_turnira = -1;
	for (int i = 0; i < TOURNAMENT_SIZE; i++) {
		int k = std::rand() % POPULATION_SIZE;
		if (populacija.at(k).fitnes > max) {
			max = populacija.at(k).fitnes;
			pobednik_turnira = k;
		}
	}
	return pobednik_turnira;
}

// funkcija koja vrsi jednopoziciono ukrstanje
void ukrstanje(Jedinka& roditelj1, Jedinka& roditelj2,
				Jedinka& dete1, Jedinka& dete2,
				std::vector<std::tuple<double,double>> releji)
{
	int n = roditelj1.kod_jedinke.size();
	int index = std::rand() % n;

	dete1.kod_jedinke.resize(0);
	dete2.kod_jedinke.resize(0);
	for (int i = 0; i < index; i++) {
		dete1.kod_jedinke.push_back(roditelj1.kod_jedinke.at(i));
		dete2.kod_jedinke.push_back(roditelj2.kod_jedinke.at(i));
	}
	for (int i = index; i < n; i++) {
		dete1.kod_jedinke.push_back(roditelj2.kod_jedinke.at(i));
		dete2.kod_jedinke.push_back(roditelj1.kod_jedinke.at(i));
	}
	dete1.popraviKodJedinke(releji);
	dete2.popraviKodJedinke(releji);
}

// funkcija koja vrsi mutaciju
void mutacija(Jedinka& jedinka,
				std::vector<std::tuple<double,double>> releji)
{
	for (int i = 0; i < jedinka.kod_jedinke.size(); i++) {
		if ((std::rand() % 1000) > MUTATION_PROB)
				jedinka.kod_jedinke.at(i) = !jedinka.kod_jedinke.at(i);
	}
	jedinka.popraviKodJedinke(releji);
}

//funkcija koja vrsi poredjenje i po njoj se jedinke sortiraju
bool poredjenje(Jedinka j1, Jedinka j2)
{
	return (j1.fitnes > j2.fitnes);
}

// funkcija koja kreira listu kordinata lokacija(tuple)
std::vector<std::tuple<int,int>> set_lokacije()
{
  std::vector<std::tuple<int,int>> lokacije = {};

  for(int i = 0; i < 287; i += 5){
    for(int j = 0; j < 287; j += 5){
      lokacije.push_back(std::tuple<int,int>(i, j));
    }
  }
  return lokacije;
}

// funkcija koja kreira listu kordinata releja(tuple)
std::vector<std::tuple<double,double>> set_releji()
{
	std::vector<std::tuple<double,double>> releji = {};

	// postavljanje 1 u centar i 2 random predajnika u okviru celije
    for(int i = 0; i < 287; i += 41){
    	for(int j = 0; j < 287; j += 41){
   			releji.push_back(std::tuple<double,double>((i + 20.5), (j + 20.5)));
              	for(int z = 0; z < 2; z++){
                    double x = (rand() % 41 +i) * 1.0;
                    double y = (rand() % 41 + j) * 1.0;
                    releji.push_back(std::tuple<double,double>(x, y));
                }
   		}
  	}
    // random mesa elemente u vektoru
	std::random_shuffle(releji.begin(), releji.end());
	return releji;
}

// funkcija koja racuna rastojanje izmedju jednog releja i jedne lokacije i vraca 1/0 u zavisnosti da li lokacija pripada tom releju
bool rastojanje(std::tuple<double,double> releji, std::tuple<int,int> lokacije)
{
	if(((std::get<0>(releji) - OPSEG) <= std::get<0>(lokacije))
			&& ((std::get<0>(releji) + OPSEG) > std::get<0>(lokacije))
			&& ((std::get<1>(releji) - OPSEG) <= std::get<1>(lokacije))
			&& ((std::get<1>(releji) + OPSEG) > std::get<1>(lokacije)))
		return true;
  	else
		return false;
}

// funkcija koja vraca mapu gde je kljuc relej a vrednost lista lokacija koje ona pokriva
std::map<std::tuple<double, double>, std::vector<std::tuple<int,int>>>
	uparivanjeLokacijaSaRelejima(std::vector<std::tuple<double,double>> releji,
								std::vector<std::tuple<int,int>> lokacije)
{
  std::map<std::tuple<double,double>, std::vector<std::tuple<int,int>>> mapa;

	for(int i = 0; i< releji.size(); i++){
		std::vector<std::tuple<int,int>> vrednost = {};
    	for(int j=0; j< lokacije.size(); j++){
    		if(rastojanje(releji.at(i), lokacije.at(j)))
				vrednost.push_back(std::tuple<int,int>(lokacije.at(j)));
    	}
		mapa.insert({ releji.at(i), vrednost });
  	}
  return mapa;
}

// funkcija koja vraca ukupan broj svih pokrivenih lokacija
int pokrivenostSvihReleja(std::map<std::tuple<double,double>,
								std::vector<std::tuple<int,int>>> mapa_releja)
{
  std::vector<std::tuple<int, int>> tmp = {};
  std::set<std::tuple<int,int>> skup;

	for(auto& i: mapa_releja){
    	std::vector<std::tuple<int, int>> tmp2 = i.second;
        tmp.insert(std::end(tmp), std::begin(tmp2), std::end(tmp2));
    }
    for(auto &elem :tmp){
   		skup.insert(elem);
    }
    return skup.size();
}

// algoritam: Standard Genetic Algorithm
void SGA(int pokrivenost, std::map<std::tuple<double,double>,
						std::vector<std::tuple<int,int>>> mapa_releja,
						std::vector<std::tuple<double,double>> releji)
{
	for (int generacija = 0; generacija < MAX_ITERATIONS; generacija++) {
		//stampanje svake iteracije
		std::sort(populacija.begin(), populacija.end(), poredjenje);
		if(generacija%1 == 0){
		std::cout << "Iteracija: " << generacija << std::endl;
		std::cout << "Fitnes: " << populacija.at(0).fitnes << std::endl;
		std::cout << "Releji: " << populacija.at(0).korisceni_releji.size() <<std::endl;
		std::cout << "Pokrivene lokacije: " << populacija.at(0).mapa_lokacija.size() << std::endl;

		}
		//elitizam
		for (int i = 0; i < ELITE_SIZE; i++) {
			nova_populacija.at(i) = populacija.at(i);
		}
		//selekcija,ukrstanje,mutacija,fitnes
		for (int i = ELITE_SIZE; i < POPULATION_SIZE; i += 2) {
			int i1 = selekcija();
			int i2 = selekcija();
			ukrstanje(populacija.at(i1), populacija.at(i2),
						nova_populacija.at(i), nova_populacija.at(i+1), releji);
			mutacija(nova_populacija.at(i), releji);
			mutacija(nova_populacija.at(i+1), releji);
			nova_populacija.at(i).fitnes =
					nova_populacija.at(i).fitnesFunkcija(pokrivenost, mapa_releja, releji);
			nova_populacija.at(i+1).fitnes =
					nova_populacija.at(i+1).fitnesFunkcija(pokrivenost, mapa_releja, releji);
		}
		populacija = nova_populacija;
	}

}



int main(int argc, char** argv)
{
    // postavljanje generatora seed-a
	srand(time(NULL));
	// definisanje vektora lokacija
	std::vector<std::tuple<int, int>> lokacije = set_lokacije();
	// definisanje vektora releja
	std::vector<std::tuple<double, double>> releji = set_releji();

	//definisanje mape releja
	std::map<std::tuple<double, double>,
			std::vector<std::tuple<int,int>>> mapa_releja
								= uparivanjeLokacijaSaRelejima(releji, lokacije);

	//racunanje pokrivenosti svih releja
	int pokrivenost = pokrivenostSvihReleja(mapa_releja);

	// kreiranje populacija za algoritam
	for (int i = 0; i < POPULATION_SIZE; i++) {
		populacija.push_back(Jedinka(mapa_releja, releji, pokrivenost));
		nova_populacija.push_back(Jedinka(mapa_releja, releji, pokrivenost));
	}

	double start = omp_get_wtime();
	//pocetak algoritma
	SGA(pokrivenost, mapa_releja, releji);
	double finish = omp_get_wtime();

	std::cout << "Vreme izvrsavanja: " << (finish - start) << std::endl;


    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);


    // Kreira se prozor.
    glutInitWindowSize(500, 500);
    glutCreateWindow(argv[0]);

     // Poziva se funkcija za inicijalizaciju
    init();

    // Registruju se callback funkcije.
    glutReshapeFunc(on_reshape);
    glutKeyboardFunc(on_keyboard);
    glutDisplayFunc(on_display);

    // Program ulazi u glavnu petlju.
    glutMainLoop();



	return 0;
}



void on_display(void)
{
    // Cistimo kolor bafer
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_QUADS);
        glColor4f(0.0, 0.0, 0.6, 0.75);
        glVertex3f(-1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, -1.0, 0.0);
        glVertex3f(-1.0, -1.0, 0.0);
    glEnd();


    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // prolazimo kroz vektor koriscenih releja i iscrtavamo
    auto najboljaJedinka = populacija.at(0).korisceni_releji;

    for(int i = 0; i<najboljaJedinka.size(); i++){
        double x = std::get<0>(najboljaJedinka.at(i));
        double y = std::get<1>(najboljaJedinka.at(i));

        // iscrtavanje opsega releja cije su koordinate x,y
        glBegin(GL_QUADS);
            glColor4f(0.6, 0.9, 0.4, 0.75);
            glVertex3f((x-OPSEG-143.5)/143.5, (y+OPSEG-143.5)/143.5, 0.0);
            glVertex3f((x+OPSEG-143.5)/143.5, (y+OPSEG-143.5)/143.5, 0.0);
            glVertex3f((x+OPSEG-143.5)/143.5, (y-OPSEG-143.5)/143.5, 0.0);
            glVertex3f((x-OPSEG-143.5)/143.5, (y-OPSEG-143.5)/143.5, 0.0);
        glEnd();

    }

    glutSwapBuffers();
}

void on_reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void on_keyboard(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:
            // Izlazimo iz programa
            exit(0);
            break;
        default:
            break;
    }

    glutPostRedisplay();
}
