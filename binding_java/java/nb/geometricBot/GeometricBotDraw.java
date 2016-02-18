/**
 * Drawing functions of the Geometric Bot
 * 2011-2015 Victor Eduardo Cardoso Nungaray
 * Twitter: @victore_cardoso
 * email: victorc@cimat.mx
 *
 */
package nb.geometricBot;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class GeometricBotDraw {
    private GeometricBotDraw() {}

    public static void drawMesh(Mesh mesh, String imageFile, int width, int height)
    {
	try {
	    /* TYPE_INT_ARGB specifies the image format: 
	     * 8-bit RGBA packed into integer pixels 
	     */
	    BufferedImage bi = 
		new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

	    Graphics2D ig2 = bi.createGraphics();

	    Box box = GeometricBot.getEnvelopingBox(mesh.vertices);
	    
	    float[] center = new float[2];
	    float zoom = GeometricBot.getCenterAndZoom(width, height, box, center);

	    /* Draw background */
	    ig2.setColor(Color.WHITE);
	    ig2.fillRect(0, 0, width, height);
	    
	    ig2.setColor(Color.BLACK);
	    drawEdges(ig2, mesh.connEdges, mesh.vertices,
		      width, height, center, zoom);
	    
	    ImageIO.write(bi, "PNG", new File(imageFile));
	    /* You could use "JPEG", "gif" and "BMP" instead of "PNG" */      
	} catch (IOException ie) {
	    ie.printStackTrace();
	}
    }

    private static void drawEdges(Graphics2D g2, int connEdges[], float vertices[],
				  int width, int height, float center[], 
				  float zoom)
    {
	for (int i = 0; i < connEdges.length / 2; i++) {
	    int id1 = connEdges[i * 2];
	    int id2 = connEdges[i*2+1];
	    float x1 =  zoom * vertices[id1 * 2] - zoom * center[0] + width / 2.0f;
	    float x2 =  zoom * vertices[id2 * 2] - zoom * center[0] + width / 2.0f;
	    float y1 =  zoom * vertices[id1*2+1] - zoom * center[1] + height / 2.0f;
	    float y2 =  zoom * vertices[id2*2+1] - zoom * center[1] + height / 2.0f;
	    g2.draw(new Line2D.Float(x1, height - y1, x2, height - y2));
	}
    }
}
