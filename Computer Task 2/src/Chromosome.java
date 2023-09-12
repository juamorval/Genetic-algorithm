import java.util.ArrayList;
import java.util.List;

public class Chromosome {
    //Chromosome =>  [2,3,0,6,4,2,7,1]  chromosome[i]=rows, i=column number
    public List<Integer> chromosome;

    public Chromosome(List<Integer> chromosome) {
        this.chromosome=chromosome;
    }

    public List<Integer> getChromosome() {
        return chromosome;
    }

    public void setChromosome(List<Integer> chromosome) {
        this.chromosome = chromosome;
    }


}
