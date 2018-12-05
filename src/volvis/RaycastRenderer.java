/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;
import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;

import java.awt.image.BufferedImage;

//import org.jetbrains.annotations.Contract;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    public static int RenderType = 0;
    private Volume volume = null;
    private GradientVolume gradients = null;
    public static boolean illumination = false;
    private double stepdefault = 1;
    private double stepfast = 10;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    //Illumination variables
    double i_a = 0.0; // ambient light
    int alpha = 10;
    double k_spec = 0.2;
    double k_diff = 0.7;
    double k_ambient = 0.7;

    int definedIntensity;
    double definedRadius;
    TFColor definedColor;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset
        //tFunc.setTestFunc();


        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }


//    short getVoxel(double[] coord) {
//
//        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
//                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
//            return 0;
//        }
//
//
//        int x = (int) Math.floor(coord[0]);
//        int y = (int) Math.floor(coord[1]);
//        int z = (int) Math.floor(coord[2]);
//        return volume.getVoxel(x,y,z);
//    }

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > (volume.getDimX() - 1) || coord[1] < 0 || coord[1] > (volume.getDimY() - 1)
                || coord[2] < 0 || coord[2] > (volume.getDimZ() - 1)) {
            return 0;
        }


        int x1 = (int) Math.floor(coord[0]);
        int y1 = (int) Math.floor(coord[1]);
        int z1 = (int) Math.floor(coord[2]);
        int x2 = (int) Math.ceil(coord[0]);
        int y2 = (int) Math.ceil(coord[1]);
        int z2 = (int) Math.ceil(coord[2]);
        double alpha = (coord[0] - x1);
        double beta = (coord[1] - y1);
        double gamma = (coord[2] - z1);
        return (short) (volume.getVoxel(x1, y1, z1) * (1 - alpha) * (1 - beta) * (1 - gamma)
                + volume.getVoxel(x2, y1, z1) * alpha * (1 - beta) * (1 - gamma)
                + volume.getVoxel(x1, y1, z2) * (1 - alpha) * (1 - beta) * gamma
                + volume.getVoxel(x2, y1, z2) * alpha * (1 - beta) * gamma
                + volume.getVoxel(x1, y2, z1) * (1 - alpha) * beta * (1 - gamma)
                + volume.getVoxel(x2, y2, z1) * alpha * beta * (1 - gamma)
                + volume.getVoxel(x1, y2, z2) * (1 - alpha) * beta * gamma
                + volume.getVoxel(x2, y2, z2) * alpha * beta * gamma);
    }

    void slicer(double[] viewMatrix) {

        // clear image
        clear_image();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();


        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);


                // BufferedImage expects a pixel color packed as ARGB in an int
                image.setRGB(i, j, toARGB(voxelColor));
            }
        }

    }

    void mip(double[] viewMatrix) {

        double step = stepdefault;
        if (interactiveMode) {
            step = stepfast;
        }
        // clear image
        clear_image();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] coeffs = new double[6];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        int maxDimension = volume.getDiagonal();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

                int val = -Integer.MAX_VALUE;
                double[] coords = new double[3];
                //for (double pos = max_neg; pos < min_pos; pos += (min_pos-max_neg)*step){
                for (double pos = -0.5 * maxDimension; pos <= 0.5 * maxDimension; pos += step) {
                    coords[0] = pixelCoord[0] + pos * viewVec[0];
                    coords[1] = pixelCoord[1] + pos * viewVec[1];
                    coords[2] = pixelCoord[2] + pos * viewVec[2];
                    val = Math.max(val, getVoxel(coords));
                }


                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                image.setRGB(i, j, toARGB(voxelColor));
            }
        }

    }

    void compositing(double[] viewMatrix) {

        double step = stepdefault;
        if (interactiveMode) {
            step = stepfast;
        }
        // clear image
        clear_image();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] coeffs = new double[6];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data

        TFColor voxelColor = new TFColor();
        TFColor previousColor = new TFColor();
        TFColor nextColor = new TFColor();

        int maxDimension = volume.getDiagonal();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
//                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
//                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
//                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
//                coeffs[0] = -pixelCoord[0]/viewVec[0];
//                coeffs[1] = -pixelCoord[1]/viewVec[1];
//                coeffs[2] = -pixelCoord[2]/viewVec[2];
//                coeffs[3] = (volume.getDimX()-pixelCoord[0])/viewVec[0];
//                coeffs[4] = (volume.getDimY()-pixelCoord[1])/viewVec[1];
//                coeffs[5] = (volume.getDimZ()-pixelCoord[2])/viewVec[2];
//                double min_pos=Double.MAX_VALUE;
//                double max_neg = -Double.MAX_VALUE;
//                for (int t = 0; t<6; ++t) {
//                    if (coeffs[t] >= 0) {
//                        min_pos = Math.min(min_pos, coeffs[t]);
//                    } else {
//                        max_neg = Math.max(max_neg, coeffs[t]);
//                    }
//                }

                int val;
                double[] coords = new double[3];
                previousColor.r = 0;
                previousColor.g = 0;
                previousColor.b = 0;
                nextColor.r = 0;
                nextColor.g = 0;
                nextColor.b = 0;
//                for (double pos = max_neg; pos < min_pos; pos += step){
                for (double pos = -0.5 * maxDimension; pos <= 0.5 * maxDimension; pos += step) {
                    coords[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + pos * viewVec[0];
                    coords[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + pos * viewVec[1];
                    coords[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + pos * viewVec[2];
                    val = getVoxel(coords);

                    voxelColor = tFunc.getColor(val);
                    if (((coords[0] < volume.getDimX() && coords[0] >= 0) && (coords[1] < volume.getDimY() && coords[1] >= 0) && (coords[2] < volume.getDimZ() && coords[2] >= 0))) {


                        nextColor.r = voxelColor.a * voxelColor.r + (1. - voxelColor.a) * previousColor.r;
                        nextColor.g = voxelColor.a * voxelColor.g + (1. - voxelColor.a) * previousColor.g;
                        nextColor.b = voxelColor.a * voxelColor.b + (1. - voxelColor.a) * previousColor.b;
                        nextColor.a = voxelColor.a + (1 - voxelColor.a) * previousColor.a;
                        //nextColor.a = 1.;

                        previousColor = nextColor;
                    }
                }
                image.setRGB(i, j, toARGB(nextColor));
            }
        }

    }

    void transfer2D(double[] viewMatrix) {

        double step = stepdefault;
        if (interactiveMode) {
            step = stepfast;
        }
        // clear image
        clear_image();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] coeffs = new double[6];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        //Colours
        TFColor voxelColor = new TFColor();
        TFColor previousColor = new TFColor();
        TFColor nextColor = new TFColor();

        // Illumination variables
        double[] lVec = new double[]{viewVec[0], viewVec[1], viewVec[2]};
        double[] hVec = new double[3]; // H-Vector
        double[] nVec = new double[3]; // N-Vector
        hVec = lVec; // Therefore L Vector is equal to H Vector
        // Assuming White color light

        // Get the GUI user defined variables
        definedIntensity = tfEditor2D.triangleWidget.baseIntensity;
        definedRadius = tfEditor2D.triangleWidget.radius;
        definedColor = tfEditor2D.triangleWidget.color;

        int maxDimension = volume.getDiagonal();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                int val;
                double[] coords = new double[3];
                previousColor.r = 1;
                previousColor.g = 1;
                previousColor.b = 1;
                nextColor.r = 0;
                nextColor.g = 0;
                nextColor.b = 0;

                for (double pos = 0.5 * maxDimension; pos >= -0.5 * maxDimension; pos -= step) {
                    coords[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + pos * viewVec[0];
                    coords[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + pos * viewVec[1];
                    coords[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + pos * viewVec[2];
                    if (((coords[0] < volume.getDimX() - 1 && coords[0] > 0)
                            && (coords[1] < volume.getDimY() - 1 && coords[1] > 0)
                            && (coords[2] < volume.getDimZ() - 1 && coords[2] >= 0))) {

                        val = getVoxel(coords);

                        VoxelGradient voxGrad = gradients.getGradient((int) Math.floor(coords[0]),
                                                                      (int) Math.floor(coords[1]),
                                                                      (int) Math.floor(coords[2]));

                        if (val == definedIntensity && voxGrad.mag == 0) {
                            voxelColor.a = definedColor.a * 1.0;
                        } else if (voxGrad.mag > tfEditor2D.triangleWidget.minMag
                                && voxGrad.mag < tfEditor2D.triangleWidget.maxMag
                                && ((val - definedRadius * voxGrad.mag) <= definedIntensity)
                                && ((val + definedRadius * voxGrad.mag) >= definedIntensity)) {
                            voxelColor.a = definedColor.a * (1.0 - (1.0 / definedRadius)
                                    * (Math.abs((definedIntensity - val) / voxGrad.mag)));
                        } else
                            voxelColor.a = 0.0;

                        if (illumination) {
                            if (voxGrad.mag >= tfEditor2D.triangleWidget.minMag
                                    && voxGrad.mag <= tfEditor2D.triangleWidget.maxMag
                                    && voxelColor.a > 0.0) {
                                // Filling N-Vector:
                                nVec[0] = voxGrad.x / voxGrad.mag;
                                nVec[1] = voxGrad.y / voxGrad.mag;
                                nVec[2] = voxGrad.z / voxGrad.mag;

                                // Computing required dot products
                                double l_dot_n = VectorMath.dotproduct(lVec, nVec);
                                double n_dot_h = VectorMath.dotproduct(nVec, hVec);

                                if (l_dot_n > 0 && n_dot_h > 0) {
                                    nextColor.r = i_a + (definedColor.r * k_diff * l_dot_n)
                                            + k_spec * Math.pow(n_dot_h, alpha);
                                    nextColor.g = i_a + (definedColor.g * k_diff * l_dot_n)
                                            + k_spec * Math.pow(n_dot_h, alpha);
                                    nextColor.b = i_a + (definedColor.b * k_diff * l_dot_n)
                                            + k_spec * Math.pow(n_dot_h, alpha);
                                    nextColor.a = 1.;
                                    break;
                                }

                            }
                        } else {

                            nextColor.r = voxelColor.a * definedColor.r + (1. - voxelColor.a) * previousColor.r;
                            nextColor.g = voxelColor.a * definedColor.g + (1. - voxelColor.a) * previousColor.g;
                            nextColor.b = voxelColor.a * definedColor.b + (1. - voxelColor.a) * previousColor.b;
                            nextColor.a = voxelColor.a + (1 - voxelColor.a) * previousColor.a;
                        }

                        previousColor = nextColor;
                    }
                }
                image.setRGB(i, j, toARGB(nextColor));
            }
        }
    }

//    private TFColor simplePhong(TFColor color_in, double[] c, boolean light, double[] viewVec) {
//            TFColor out = new TFColor(color_in.r, color_in.g, color_in.b, color_in.a);
//
//            //Calculate normal according to gradients
//            double[] N = new double[3];
//            VoxelGradient gradient = gradients.getGradient((int) Math.floor(c[0]), (int) Math.floor(c[1]), (int) Math.floor(c[2]));
//
//            if (gradient.mag > 0.0 && color_in.a > 0.0) {
//                // Set L (= V) = H to be the vector pointing from the point to our 'eye'/light source
//                double[] lVec = new double[]{viewVec[0], viewVec[1], viewVec[2]};
//                double[] hVec = new double[3]; // H-Vector
//                double[] nVec = new double[3]; // N-Vector
//
//                // Calculate normal vector N
//                nVec[0] = (double) gradient.x / (double) gradient.mag;
//                nVec[1] = (double) gradient.y / (double) gradient.mag;
//                nVec[2] = (double) gradient.z / (double) gradient.mag;
//
//                // Compute L dot N and N dot H
//                double ln = VectorMath.dotproduct(lVec, nVec);
//                double nh = VectorMath.dotproduct(nVec, hVec);
//
//                if (ln > 0 && nh > 0) {
//                    out.r = i_a + color_in.r * (k_diff * ln) + k_spec * Math.pow(nh, alpha);
//                    out.g = i_a + color_in.g * (k_diff * ln) + k_spec * Math.pow(nh, alpha);
//                    out.b = i_a + color_in.b * (k_diff * ln) + k_spec * Math.pow(nh, alpha);
//                }
//            }
//            return out;
//    }

    private void clear_image() {
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }

    private int toARGB(TFColor in) {
        int c_alpha = in.a <= 1.0 ? (int) Math.floor(in.a * 255) : 255;
        int c_red = in.r <= 1.0 ? (int) Math.floor(in.r * 255) : 255;
        int c_green = in.g <= 1.0 ? (int) Math.floor(in.g * 255) : 255;
        int c_blue = in.b <= 1.0 ? (int) Math.floor(in.b * 255) : 255;
        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        return pixelColor;
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        switch (RenderType) {
            case 0:
                slicer(viewMatrix);
                break;
            case 1:
                mip(viewMatrix);
                break;
            case 2:
                compositing(viewMatrix);
                break;
            case 3:
                transfer2D(viewMatrix);
                break;
            default:
                System.out.println("Invalid renderer method.");
                break;
        }
        if (RenderType == 0)
            slicer(viewMatrix);
        else if (RenderType == 1)
            mip(viewMatrix);
        else if (RenderType == 2)
            compositing(viewMatrix);

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }

    // Method for getting the maximum value
    private static int getMax(int[] inputArray) {

        int maxValue = inputArray[0];

        for (int i = 1; i < inputArray.length; i++) {

            if (inputArray[i] > maxValue) {

                maxValue = inputArray[i];
            }
        }

        return maxValue;
    }

    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
