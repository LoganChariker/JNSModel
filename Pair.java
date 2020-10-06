/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package sim;

/**
 *
 * @author charikar
 */
public class Pair implements Comparable<Pair> {
    public int x;
    public int index;
    
    public Pair(int x, int index){
        this.x = x;
        this.index = index;
    }   

    public int compareTo(Pair p2) {
        if (this.x > p2.x){
            return 1;
        } else if (this.x < p2.x){
            return -1;
        } else {
            return 0;
        }
    }

}
