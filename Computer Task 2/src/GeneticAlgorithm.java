import java.util.*;

public class GeneticAlgorithm {
    public static int Nqueens = 20;             //queen's number
    public static int popSize = 200;            //population size
    //public int stop_cicles;                  //stop condition
    public static int maxClashes = Nqueens * (Nqueens-1)/2;
    public static double elitism = maxClashes - Nqueens;       //elitism probability
    public static double pc = 3.0;             //crossover probability
    public static double pm = 3.0;             //mutation probability
    public static int offspring = 60;          //n offspring created
    public static int chosen = 10;              //chosen subset size in tournament selection

    /*
    public List<Chromosome> population;
    public static Map<Double, Chromosome> selected;
    public List<Chromosome> new_population;

     */

    public static Chromosome algorithm(Integer popSize, Integer n, Integer nOffspring, Double elitism,
                                 Integer chosen, Double pc, Double pm, Integer maxClashes) {

        List<Chromosome> population = initPopulation(popSize, n);
        int c = 1;
        System.out.println("Evolutionary cycle: 0");
        List<List<Integer>> generation = new ArrayList<>();
        for (int i=0; i<population.size(); i++) {
            generation.add(population.get(i).getChromosome());

        }

        List<Chromosome> bests = bestIndividuals(population, elitism);

        List<Chromosome> newPopulation = new ArrayList<>();
        population.clear();
        for (int i = 0; i < bests.size(); i++) {
            population.add(bests.get(i));
        }


        System.out.println("Population:");
        System.out.println(generation);
        Chromosome winner = null;

        boolean flag = true;
        while(flag) {

            for(int h=newPopulation.size(); h<nOffspring; h=newPopulation.size()) {
                Chromosome parentA = tournamentSelection(population, chosen);
                Chromosome parentB = tournamentSelection(population, chosen);
                newPopulation.add(parentA);
                newPopulation.add(parentB);

                List<Chromosome> children = crossover(parentA, parentB, pc);

                Chromosome childA = mutation(children.get(0), pm);
                Chromosome childB = mutation(children.get(1), pm);

                if(fitness(childA) > fitness(childB)) {
                    newPopulation.add(childA);
                } else {
                    newPopulation.add(childB);
                }

            }



            //Replacement
            population.clear();
            for(int i=0; i< newPopulation.size(); i++) {
                population.add(newPopulation.get(i));
            }
            newPopulation.clear();


            System.out.println("\nEvolutionary cycle: " + c);

            generation.clear();
            for(int j=0; j< population.size(); j++) {
                generation.add(population.get(j).getChromosome());
            }
            System.out.println("Population: ");
            System.out.println(generation);

            int i=0;
            while (i < population.size()) {
                if (maxClashes==fitness(population.get(i))) {
                    winner = population.get(i);
                    flag = false;
                    break;
                }
                i++;
            }

            c++;

        }

        return winner;

    }



    public static Chromosome createChromosome(Integer n) {
        List<Integer> res = new ArrayList<>();
        for(int i=0; i<n; i++) {
            res.add((int) (Math.random() * n));
        }

        Chromosome chr = new Chromosome(res);

        return chr;
    }

    public static List<Chromosome> initPopulation(Integer popSize, Integer n) {
        List<Chromosome> population = new ArrayList<>();
        for(int i=0; i<popSize; i++) {
            Chromosome chromosome = createChromosome(n);
            population.add(chromosome);
        }

        return population;
    }

    public static double fitness(Chromosome chromosome) {  //chromosome defined to do not have vertical clashes
        double res = 0;
        List<Integer> chr = chromosome.getChromosome();
        for(int i=0; i<chr.size(); i++) {  //i->row of i
            for(int j=i+1; j<chr.size(); j++) {  //j->comparing to the other rows
                if(chr.get(i).equals(chr.get(j)) ||  //horizontal clash
                        Math.abs(chr.get(i)-chr.get(j))==Math.abs(i-j)) {   //diagonal clash
                    res++;
                }
            }
        }
        double maxClash = chr.size() * (chr.size()-1)/2;     // N=8; maxClash=7+6+5+4+3+2+1 (all queens in the same row)
        double fitness = maxClash-res;                      //maxClash - nº collisions
        return fitness;                                     //better fitness = better chromosomes

    }

    public static List<Chromosome> bestIndividuals(List<Chromosome> population, Double elitism) {     //piojo's eye
        List<Chromosome> selected = new ArrayList<>();

        for(int i=0; i<population.size(); i++) {
            if((fitness(population.get(i)))>elitism) {
                selected.add(population.get(i));
            }
        }

        return selected;

    }

    public static Chromosome tournamentSelection(List<Chromosome> bests, Integer chosen) {
        List<Chromosome> subset = new ArrayList<>();
        Double maxFitness = 0.0;
        Chromosome father = null;
        List<Integer> lastRandom = new ArrayList<>();

        while(chosen!= subset.size()) {
            int random = (int) (Math.random()* bests.size());
            if(!lastRandom.contains(random)) {
                subset.add(bests.get(random));
                lastRandom.add(random);
            }

        }


        for(int i=0; i<subset.size(); i++) {
            if(fitness(subset.get(i))>maxFitness) {
                father = subset.get(i);
                maxFitness=fitness(subset.get(i));
            }
        }

        return father;

    }

    public static List<Chromosome> crossover(Chromosome parentA, Chromosome parentB, Double pc) {
        List<Chromosome> res = new ArrayList<>();
        List<Integer> decodeA = parentA.getChromosome();
        List<Integer> decodeB = parentB.getChromosome();

        List<Integer> childA = new ArrayList<>();
        List<Integer> childB = new ArrayList<>();

        for(int i=0; i<pc; i++) {
            childA.add(decodeA.get(i));
            childB.add(decodeB.get(i));
        }

        for(double i=pc; i<decodeA.size(); i++) {
            childA.add(decodeB.get((int) i));
            childB.add(decodeA.get((int) i));
        }

        Chromosome chA = new Chromosome(childA);
        Chromosome chB = new Chromosome(childB);

        res.add(chA);
        res.add(chB);

        return res;

    }

    public static Chromosome mutation(Chromosome chromosome, Double pm) {
        List<Integer> decode = chromosome.getChromosome();

        for(int i=0; i<decode.size(); i++) {
            double random = Math.random()* decode.size();
            if(random<pm) {
                double mutation = Math.random()* decode.size();
                decode.set(i, (int) mutation);

            }
        }

        Chromosome mutated = new Chromosome(decode);

        return mutated;

    }

    public static void printChessboard(Chromosome chr) {
        int c = 0;
        for(int i=0; i<chr.getChromosome().size(); i++) {       //row
            for(int j=0; j<chr.getChromosome().size(); j++) {   //column
                if(chr.getChromosome().get(j)==i) {
                    System.out.print("♛ ");
                    c++;
                } else {
                    if(c%2==0) {
                        System.out.print("▢ ");
                        c++;
                    } else {
                        System.out.print("▨ ");
                        c++;
                    }

                }
            }
            System.out.println();
            c++;
        }
    }


    public static void main(String args[]) {
        Chromosome winner = algorithm(popSize, Nqueens, offspring, elitism, chosen, pc, pm, maxClashes);
        System.out.println("\nThe winner chromosome is: " +
                winner.getChromosome() + "\n");
        printChessboard(winner);
    }
}
