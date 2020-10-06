/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package sim;

/**
 *
 * @author charikar
 */
public class NetHelper {
public NetHelper(){
        System.out.println("You've created a NetHelper.  Why did you do this?");
    }

    /* convert a location in cortex to an equivalent location in the first hypercolumn with regards to orientation, making the center the origin */
    public static double movCoordToHC1(double xory, double HCWidth){
        double s = xory / HCWidth;
        double m = s % 1.0;
        int flipHuh = (((int)Math.floor(s)) % 2);
        if (flipHuh == 1){
            /* flip */
            return (1.0 - m) - 0.5;
        } else {
            return m - 0.5;
        }
    }

    public static MyPoint movLocToHC1(double x, double y, double HCWidth){
        return new MyPoint(movCoordToHC1(x,HCWidth),movCoordToHC1(y,HCWidth));
    }

    /* corresponds a cortex location to an orientation */
    public static MyPoint locToOri(MyPoint loc, double HCWidth){
        MyPoint HC1Loc = movLocToHC1(loc.x,loc.y, HCWidth);

        HC1Loc.normalize();

        double x = HC1Loc.x;
        double y = HC1Loc.y;
        if (y > 0.001 || y < -0.001) {
            MyPoint result = new MyPoint(0.5*(x+1.0), -0.5*y);
            result.normalize();
            return result;
        } else if (x > 0.0) {
            return new MyPoint(1.0,0.0);
        } else {
            return new MyPoint(0.0,1.0);
        }
    }

    /* give the angle from the positive x-axis of a point in R^2 */
    public static double pointToAngle(double x, double y){
        double r = Math.sqrt(x*x+y*y);

        double angle;
        if (y > 0){
            angle = Math.acos(x / r);
        } else if (y == 0 && x < 0){
            angle = Math.PI;
        } else {
            angle = 2.0*Math.PI-Math.acos(x / r);
        }

        return angle;
    }
}
