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
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    // set the mode view:
    // 0: slicer
    // 1: mip view
    // 2: copositing
    // 3: 2d without shading
    // 4: 2d whith shading
    //default = slicer()
    public int mode = 0;

    public double resolution = 15.0;

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

    boolean inOfBox(double[] coord) {
        // check if the coordinates are inside the bounding box
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY() || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return false;
        } else {
            return true;
        }
    }

    short getVoxel(double[] coord) {

        if (!inOfBox(coord)) {
            return 0;
        }

        //int x = (int) Math.floor(coord[0]);
        //int y = (int) Math.floor(coord[1]);
        //int z = (int) Math.floor(coord[2]);
        //return volume.getVoxel(x, y, z);
        try {

            //set the original cood as the letter without number
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];

            //set the "floor" int to the lowest cord
            int x0 = (int) Math.floor(coord[0]);
            int y0 = (int) Math.floor(coord[1]);
            int z0 = (int) Math.floor(coord[2]);

            //set the "ceil" int to the higest cord
            int x1 = (int) Math.ceil(coord[0]);
            int y1 = (int) Math.ceil(coord[1]);
            int z1 = (int) Math.ceil(coord[2]);

            //foreach axis compute the "apha" (alpha, beta and gamma) as the formula in slide 6
            double alpha = (x - x0) / (x1 - x0);
            double beta = (y - y0) / (y1 - y0);
            double gamma = (z - z0) / (z1 - z0);

            //get the eight values sourronding the voxel
            short sx0 = (short) volume.getVoxel(x1, y1, z1);
            short sx1 = (short) volume.getVoxel(x1 + 1, y1, z1);
            short sx2 = (short) volume.getVoxel(x1, y1 + 1, z1);
            short sx3 = (short) volume.getVoxel(x1, y1, z1 + 1);
            short sx4 = (short) volume.getVoxel(x1 + 1, y1, z1 + 1);
            short sx5 = (short) volume.getVoxel(x1, y1 + 1, z1 + 1);
            short sx6 = (short) volume.getVoxel(x1 + 1, y1 + 1, z1 + 1);
            short sx7 = (short) volume.getVoxel(x1 + 1, y1 + 1, z1);

            // computhe sx (as formula in slide 7) and return the value
            double sx = (double) ((1 - alpha) * (1 - beta) * (1 - gamma) * sx0) + ((alpha) * (1 - beta) * (1 - gamma) * sx1) + ((alpha) * (beta) * (1 - gamma) * sx2) + ((1 - alpha) * (beta) * (1 - gamma) * sx3) + ((1 - alpha) * (1 - beta) * (gamma) * sx4) + ((alpha) * (1 - beta) * (gamma) * sx5) + ((alpha) * (beta) * (gamma) * sx6) + ((1 - alpha) * (beta) * (gamma) * sx7);
            return (short) sx;

        } catch (Exception e) {
            // if it is out of bounds, just return 0
            return 0;
        }

    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

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

                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

                int val = getVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

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

        // calculate diagonal 
        int vectorLenght = (int) (Math.sqrt(Math.pow(volume.getDimX(), 2) + Math.pow(volume.getDimY(), 2) + Math.pow(volume.getDimZ(), 2)));
        int vectorCenter = (int) vectorLenght / 2;

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                int maximum = 0;

                double amount;
                // if image is static, use full resolution
                if (!interactiveMode) {
                    amount = 1.0;
                } else {
                    // if image is in movement, use full less resolution
                    amount = resolution;
                }

                // third loop for the depth of the vector
                for (double k = -vectorCenter; k < vectorCenter; k = k + amount) {

                    //pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0] + viewVec[0] * k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1] + viewVec[1] * k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2] + viewVec[2] * k;

                    int val = getVoxel(pixelCoord);
                    //if maximum is not maximum, set it
                    if (val > maximum) {
                        maximum = val;
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maximum / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maximum > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);

            }
        }
    }

    void compositing(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

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

        TFColor voxelColor = new TFColor();

        // calculate diagonal 
        int vectorLenght = (int) (Math.sqrt(Math.pow(volume.getDimX(), 2) + Math.pow(volume.getDimY(), 2) + Math.pow(volume.getDimZ(), 2)));
        int vectorCenter = (int) vectorLenght / 2;

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                int val = 0;
                TFColor prevColor = new TFColor(0, 0, 0, 1);

                double amount;
                // if image is static, use full resolution
                if (!interactiveMode) {
                    amount = 1.0;
                } else {
                    // if image is in movement, use full less resolution
                    amount = resolution;
                }

                for (double k = -vectorCenter; k < vectorCenter; k = k + amount) {

                    //pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0] + viewVec[0] * k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1] + viewVec[1] * k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2] + viewVec[2] * k;

                    val = getVoxel(pixelCoord);

                    TFColor userColor = tFunc.getColor(val);

                    // set new colors using back-to-front compositing order (slide 3-3)
                    voxelColor.r = (userColor.r * userColor.a) + ((1 - userColor.a) * voxelColor.r);
                    voxelColor.g = (userColor.g * userColor.a) + ((1 - userColor.a) * voxelColor.g);
                    voxelColor.b = (userColor.b * userColor.a) + ((1 - userColor.a) * voxelColor.b);
                    voxelColor.a = (1 - userColor.a) * voxelColor.a;

                    prevColor.r = voxelColor.r;
                    prevColor.g = voxelColor.g;
                    prevColor.b = voxelColor.b;
                    prevColor.a = voxelColor.a;
                }

                
                // levoy's formula for alpha
                double alpha = 1 - (prevColor.a * (1 - voxelColor.a));
                
                // Map the intensity to a grey value by linear scaling
                //voxelColor.r = maximum/max;
                //voxelColor.g = voxelColor.r;
                //voxelColor.b = voxelColor.r;
                //voxelColor.a = maximum > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = alpha <= 1.0 ? (int) Math.floor(alpha * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);

            }
        }
    }

    void twodimension(double[] viewMatrix, boolean illuminated) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

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
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                // create new colors to handle operations
                TFColor currentColor = new TFColor(0, 0, 0, 1);
                TFColor prevColor = new TFColor(0, 0, 0, 1);

                // calculate diagonal 
                int vectorLenght = (int) (Math.sqrt(Math.pow(volume.getDimX(), 2) + Math.pow(volume.getDimY(), 2) + Math.pow(volume.getDimZ(), 2)));
                int vectorCenter = (int) vectorLenght / 2;

                double amount;
                // if image is static, use full resolution
                if (!interactiveMode) {
                    amount = 1.0;
                } else {
                    // if image is in movement, use full less resolution
                    amount = resolution;
                }

                for (double k = -vectorCenter; k < vectorCenter; k = k + amount) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0] + viewVec[0] * k;
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1] + viewVec[1] * k;
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2] + viewVec[2] * k;

                    int val = getVoxel(pixelCoord);

                    //Avoid errors by checking if it is inside the bounding box
                    // and if value is different from 0 to avoid null pointer exception
                    if ((inOfBox(pixelCoord)) && val != 0) {

                        //get gradient parameters
                        float magnitude = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]).mag;
                        float gradientX = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]).x;
                        float gradientY = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]).y;
                        float gradientZ = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]).z;
                        double radius = tfEditor2D.triangleWidget.radius;
                        int baseIntensity = tfEditor2D.triangleWidget.baseIntensity;

                        // get the color predefined by the user
                        TFColor userColor = tfEditor2D.triangleWidget.color;

                        // Compute opacity using the levoy's formula
                        if (magnitude == 0 && val == baseIntensity) {
                            voxelColor.a = userColor.a;
                        } else if ((magnitude > 0) && (val - radius * Math.abs(magnitude) <= baseIntensity) && (val + radius * Math.abs(magnitude) >= baseIntensity)) {
                            voxelColor.a = userColor.a * (1 - (1 / radius) * (Math.abs((baseIntensity - val) / Math.abs(magnitude))));
                        } else {
                            voxelColor.a = 0;
                        }

                        // Check if the ilumination flag is activated.
                        // if the user wants ilumintion, apply the Phong functions
                        // also check alpha and magnitude are positive to avoid errrors
                        if ((illuminated) && (voxelColor.a > 0) && (magnitude > 0)) {

                            // default parameters for the phong function
                            double kAmbient = 0.1;
                            double kDiff = 0.7;
                            double kSpec = 0.2;
                            int phongAlpha = 10;

                            // Vector of observation
                            // We need to multiply by -1 because light must go in direction to observer
                            double[] V = {-viewVec[0], -viewVec[1], -viewVec[2]};

                            // Vector in direction of maximum Highlight (levoy's)
                            // we assume H = V
                            double[] H = new double[3];
                            VectorMath.setVector(H, V[0], V[1], V[2]);

                            // assuming we have a headlight L=V (slide 3-15), so L=V
                            double[] L = new double[3];
                            VectorMath.setVector(L, V[0], V[1], V[2]);

                            // Surface normal at the position of the voxel
                            double[] N = new double[3];
                            VectorMath.setVector(N, gradientX / magnitude, gradientY / magnitude, gradientZ / magnitude);

                            //LN = NH
                            double LN = VectorMath.dotproduct(L, N);
                            double NH = VectorMath.dotproduct(N, H);

                            // compute the colors using the standard Phong model (slide 2-40)
                            voxelColor.r = kAmbient + (userColor.r * kDiff * LN) + (kSpec * Math.pow(NH, phongAlpha));
                            voxelColor.g = kAmbient + (userColor.g * kDiff * LN) + (kSpec * Math.pow(NH, phongAlpha));
                            voxelColor.b = kAmbient + (userColor.b * kDiff * LN) + (kSpec * Math.pow(NH, phongAlpha));

                            // set new colors using back-to-front compositing order (slide 3-3)
                            currentColor.r = (voxelColor.a * voxelColor.r) + ((1 - voxelColor.a) * prevColor.r);
                            currentColor.g = (voxelColor.a * voxelColor.g) + ((1 - voxelColor.a) * prevColor.g);
                            currentColor.b = (voxelColor.a * voxelColor.b) + ((1 - voxelColor.a) * prevColor.b);
                            currentColor.a = (1 - voxelColor.a) * prevColor.a;
                            
                            // set the old color as the new color
                            prevColor.r = currentColor.r;
                            prevColor.g = currentColor.g;
                            prevColor.b = currentColor.b;
                            prevColor.a = currentColor.a;

                        } else {
                            // set new colors using back-to-front compositing order (slide 3-3)
                            currentColor.r = (voxelColor.a * userColor.r) + ((1 - voxelColor.a) * prevColor.r);
                            currentColor.g = (voxelColor.a * userColor.g) + ((1 - voxelColor.a) * prevColor.g);
                            currentColor.b = (voxelColor.a * userColor.b) + ((1 - voxelColor.a) * prevColor.b);                        
                            currentColor.a = (1 - voxelColor.a) * prevColor.a;
                            
                            // set the old color as the new color
                            prevColor.r = currentColor.r;
                            prevColor.g = currentColor.g;
                            prevColor.b = currentColor.b;
                            prevColor.a = currentColor.a;
                        }
                    }
                }

                // levoy's formula for alpha
                double alpha = 1 - (prevColor.a * (1 - voxelColor.a));

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = alpha <= 1.0 ? (int) Math.floor(alpha * 255) : 255;
                int c_red = currentColor.r <= 1.0 ? (int) Math.floor(currentColor.r * 255) : 255;
                int c_green = currentColor.g <= 1.0 ? (int) Math.floor(currentColor.g * 255) : 255;
                int c_blue = currentColor.b <= 1.0 ? (int) Math.floor(currentColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
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

        if (mode == 0) {
            slicer(viewMatrix);
        } else if (mode == 1) {
            mip(viewMatrix);
        } else if (mode == 2) {
            compositing(viewMatrix);
        } else if (mode == 3) {
            twodimension(viewMatrix, false);
        } else if (mode == 4) {
            twodimension(viewMatrix, true);
        }

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
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
