/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package sim;

import java.awt.geom.Point2D;

/**
 *
 * @author Logan
 */
public class MyPoint extends Point2D.Double {

    public MyPoint(double x, double y) {
        this.x=x;
        this.y=y;
    }

    /*
    public double distTo(MyPoint a) {
        return Math.sqrt(Math.pow(this.x - a.x, 2) + Math.pow(this.y - a.y, 2));
    }

    public double distToSqr(MyPoint a) {
        double dx = this.x - a.x;
        double dy = this.y - a.y;
        return dx * dx + dy * dy;
    }
    
    
    public double normSqr() {
        return x * x + y * y;
    }

    
    
    public MyPoint rotateBy(double theta) {
        return new MyPoint(Math.cos(theta) * x - Math.sin(theta) * y, Math.sin(theta) * x + Math.cos(theta) * y);
    }
    */
    
    public double innerProduct(MyPoint a) {
        return x * a.x + y * a.y;
    }
    
    public MyPoint translateCopy(double shiftX, double shiftY){
        return new MyPoint(x+shiftX,y+shiftY);
    }

    public double getNorm(){
        return distance(0.0,0.0);
    }

    public void normalize(){
        double norm = getNorm();
        x /= norm;
        y /= norm;
    }
    
    public MyPoint scalarMultCopy(double c){
        return new MyPoint(c*x,c*y);
    }
    
    public double corticalDistance(MyPoint a, double patchWidth){
        double ax = a.x;
        double ay = a.y;
        return Math.min(distance(ax,ay),
                Math.min(distance(ax + sgn(x-ax)*patchWidth, ay),
                        Math.min(distance(ax,ay + sgn(y-ay)*patchWidth),
                                distance(ax + sgn(x-ax)*patchWidth, ay + sgn(y-ay)*patchWidth))));
    }

    public double corticalDistanceSq(MyPoint a, double patchWidth){
        double d = corticalDistance(a,patchWidth);
        return d*d;
    }
    
    private double sgn(double x){
        if (x >= 0) return 1.0; else return -1.0;
    }
}