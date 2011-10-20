##  ImageJ/NIH Image 64 byte ROI outline header
##     2 byte numbers are big-endian signed shorts
##     0-3     "Iout"
##     4-5     version (>=217)
##     6-7     roi type
##     8-9     top
##     10-11   left
##     12-13   bottom
##     14-15   right
##     16-17   NCoordinates
##     18-33   x1,y1,x2,y2 (straight line)
##     34-35   stroke width (v1.43i or later)
##     36-39   ShapeRoi size (type must be 1 if this value>0)
##     40-43   stroke color (v1.43i or later)
##     44-47   fill color (v1.43i or later)
##     48-49   subtype (v1.43k or later)
##     50-51   options (v1.43k or later)
##     52-52   arrow style or aspect ratio (v1.43p or later)
##     53-53   arrow head size (v1.43p or later)
##     54-55   rounded rect arc size (v1.43p or later)
##     56-59   position
##     60-63   reserved (zeros)
##     64-       x-coordinates (short), followed by y-coordinates
##


getByte <- function(con) {
  return(readBin(con, raw(0), 1))
}

getShort <- function(con) {
  n <- readBin(con, integer(0), 2, signed=TRUE)
  if (n < -5000) {
    seek(con, -2, origin="current")
    n <- readBin(con, integer(0), 2, signed=FALSE)
  }
  return(n)
}
    
    int getInt(int base) {
        int b0 = data[base]&255;
        int b1 = data[base+1]&255;
        int b2 = data[base+2]&255;
        int b3 = data[base+3]&255;
        return ((b0<<24) + (b1<<16) + (b2<<8) + b3);
    }

    float getFloat(int base) {
        return Float.intBitsToFloat(getInt(base));
    }


##
read.roi <- function(path) {
  ## offsets
  VERSION_OFFSET = 4;
  TYPE = 6;
  TOP = 8;
  LEFT = 10;
  BOTTOM = 12;
  RIGHT = 14;
  N_COORDINATES = 16;
  X1 = 18;
  Y1 = 22;
  X2 = 26;
  Y2 = 30;
  STROKE_WIDTH = 34;
  SHAPE_ROI_SIZE = 36;
  STROKE_COLOR = 40;
  FILL_COLOR = 44;
  SUBTYPE = 48;
  OPTIONS = 50;
  ARROW_STYLE = 52;
  ELLIPSE_ASPECT_RATIO = 52;
  ARROW_HEAD_SIZE = 53;
  ROUNDED_RECT_ARC_SIZE = 54;
  POSITION = 56;
  COORDINATES = 64;
  
  ## subtypes
  TEXT = 1;
  ARROW = 2;
  ELLIPSE = 3;
  
  ## options
  SPLINE_FIT = 1;
  DOUBLE_HEADED = 2;
  OUTLINE = 4;
    
  ## types
  polygon=0
  rect=1
  oval=2
  line=3
  freeline=4
  polyline=5
  noRoi=6
  freehand=7
  traced=8
  angle=9
  point=10;

  ## Main code
  if (!is.null(path)) {
    size <- file.info(path)$size
    if (!grepl(".roi$", path) && size>5242880)
      stop("This is not an ROI or file size>5MB)")
    ## FIXME name = f.getName();

  }
  con <- file(path)
  ##  data <- readBin(path, raw(0), size)

  ## int total = 0;
  ## while (total<size)
  ##   total += is.read(data, total, size-total);
  ## is.close();
  
  if (getByte(con)!=73 || getByte(con) != 111) {  ## "Iout"
    stop("This is not an ImageJ ROI");
  }
  version = getShort(VERSION_OFFSET);
type = getByte(TYPE);
subtype = getShort(SUBTYPE);
top= getShort(TOP);
left = getShort(LEFT);
bottom = getShort(BOTTOM);
right = getShort(RIGHT);
width = right-left;
height = bottom-top;
n = getShort(N_COORDINATES);
options = getShort(OPTIONS);
position = getInt(POSITION);
        
        if (name!=null && name.endsWith(".roi"))
            name = name.substring(0, name.length()-4);
        boolean isComposite = getInt(SHAPE_ROI_SIZE)>0;
        
        Roi roi = null;
        if (isComposite) {
            roi = getShapeRoi();
            if (version>=218) getStrokeWidthAndColor(roi);
            roi.setPosition(position);
            return roi;
        }

        switch (type) {
            case rect:
                roi = new Roi(left, top, width, height);
                int arcSize = getShort(ROUNDED_RECT_ARC_SIZE);
                if (arcSize>0)
                    roi.setCornerDiameter(arcSize);
                break;
            case oval:
                roi = new OvalRoi(left, top, width, height);
                break;
            case line:
                int x1 = (int)getFloat(X1);     
                int y1 = (int)getFloat(Y1);     
                int x2 = (int)getFloat(X2);     
                int y2 = (int)getFloat(Y2);
                if (subtype==ARROW) {
                    roi = new Arrow(x1, y1, x2, y2);        
                    ((Arrow)roi).setDoubleHeaded((options&DOUBLE_HEADED)!=0);
                    ((Arrow)roi).setOutline((options&OUTLINE)!=0);
                    int style = getByte(ARROW_STYLE);
                    if (style>=Arrow.FILLED && style<=Arrow.OPEN)
                        ((Arrow)roi).setStyle(style);
                    int headSize = getByte(ARROW_HEAD_SIZE);
                    if (headSize>=0 && style<=30)
                        ((Arrow)roi).setHeadSize(headSize);
                } else
                    roi = new Line(x1, y1, x2, y2);     
                //IJ.write("line roi: "+x1+" "+y1+" "+x2+" "+y2);
                break;
            case polygon: case freehand: case traced: case polyline: case freeline: case angle: case point:
                    //IJ.write("type: "+type);
                    //IJ.write("n: "+n);
                    //IJ.write("rect: "+left+","+top+" "+width+" "+height);
                    if (n==0) break;
                    int[] x = new int[n];
                    int[] y = new int[n];
                    int base1 = COORDINATES;
                    int base2 = base1+2*n;
                    int xtmp, ytmp;
                    for (int i=0; i<n; i++) {
                        xtmp = getShort(base1+i*2);
                        if (xtmp<0) xtmp = 0;
                        ytmp = getShort(base2+i*2);
                        if (ytmp<0) ytmp = 0;
                        x[i] = left+xtmp;
                        y[i] = top+ytmp;
                        //IJ.write(i+" "+getShort(base1+i*2)+" "+getShort(base2+i*2));
                    }
                    if (type==point) {
                        roi = new PointRoi(x, y, n);
                        break;
                    }
                    int roiType;
                    if (type==polygon)
                        roiType = Roi.POLYGON;
                    else if (type==freehand) {
                        roiType = Roi.FREEROI;
                        if (subtype==ELLIPSE) {
                            double ex1 = getFloat(X1);      
                            double ey1 = getFloat(Y1);      
                            double ex2 = getFloat(X2);      
                            double ey2 = getFloat(Y2);
                            double aspectRatio = getFloat(ELLIPSE_ASPECT_RATIO);
                            roi = new EllipseRoi(ex1,ey1,ex2,ey2,aspectRatio);
                            break;
                        }
                    } else if (type==traced)
                        roiType = Roi.TRACED_ROI;
                    else if (type==polyline)
                        roiType = Roi.POLYLINE;
                    else if (type==freeline)
                        roiType = Roi.FREELINE;
                    else if (type==angle)
                        roiType = Roi.ANGLE;
                    else
                        roiType = Roi.FREEROI;
                    roi = new PolygonRoi(x, y, n, roiType);
                    break;
            default:
                throw new IOException("Unrecognized ROI type: "+type);
        }
        if (name!=null) roi.setName(name);
        
        // read stroke width, stroke color and fill color (1.43i or later)
        if (version>=218) {
            getStrokeWidthAndColor(roi);
            boolean splineFit = (options&SPLINE_FIT)!=0;
            if (splineFit && roi instanceof PolygonRoi)
                ((PolygonRoi)roi).fitSpline();
        }
        
        if (version>=218 && subtype==TEXT)
            roi = getTextRoi(roi);

        roi.setPosition(position);
        return roi;
    }

}


## Decodes an ImageJ, NIH Image or Scion Image ROI file. 
readroi <- function(file) {
    
    private byte[] data;
    private String path;
    private String name;
    private int size;
  con <- file(file, "rb")
  bytes <- readBin(con, 
    /** Constructs an RoiDecoder using a byte array. */
    public RoiDecoder(byte[] bytes, String name) {
        is = new ByteArrayInputStream(bytes);   
        this.name = name;
        this.size = bytes.length;
    }

    
    void getStrokeWidthAndColor(Roi roi) {
        int strokeWidth = getShort(STROKE_WIDTH);
        if (strokeWidth>0)
            roi.setStrokeWidth(strokeWidth);
        int strokeColor = getInt(STROKE_COLOR);
        if (strokeColor!=0) {
            int alpha = (strokeColor>>24)&0xff;
            roi.setStrokeColor(new Color(strokeColor, alpha!=255));
        }
        int fillColor = getInt(FILL_COLOR);
        if (fillColor!=0) {
            int alpha = (fillColor>>24)&0xff;
            roi.setFillColor(new Color(fillColor, alpha!=255));
        }
    }

    public Roi getShapeRoi() throws IOException {
        int type = getByte(TYPE);
        if (type!=rect)
            throw new IllegalArgumentException("Invalid composite ROI type");
        int top= getShort(TOP);
        int left = getShort(LEFT);
        int bottom = getShort(BOTTOM);
        int right = getShort(RIGHT);
        int width = right-left;
        int height = bottom-top;
        int n = getInt(SHAPE_ROI_SIZE);

        ShapeRoi roi = null;
        float[] shapeArray = new float[n];
        int base = COORDINATES;
        for(int i=0; i<n; i++) {
            shapeArray[i] = getFloat(base);
            base += 4;
        }
        roi = new ShapeRoi(shapeArray);
        if (name!=null) roi.setName(name);
        return roi;
    }
    
    Roi getTextRoi(Roi roi) {
        Rectangle r = roi.getBounds();
        int hdrSize = RoiEncoder.HEADER_SIZE;
        int size = getInt(hdrSize);
        int style = getInt(hdrSize+4);
        int nameLength = getInt(hdrSize+8);
        int textLength = getInt(hdrSize+12);
        char[] name = new char[nameLength];
        char[] text = new char[textLength];
        for (int i=0; i<nameLength; i++)
            name[i] = (char)getShort(hdrSize+16+i*2);
        for (int i=0; i<textLength; i++)
            text[i] = (char)getShort(hdrSize+16+nameLength*2+i*2);
        Font font = new Font(new String(name), style, size);
        Roi roi2 = new TextRoi(r.x, r.y, new String(text), font);
        roi2.setStrokeColor(roi.getStrokeColor());
        roi2.setFillColor(roi.getFillColor());
        return roi2;
    }

    
    /** Opens an ROI from a byte array. */
    public static Roi openFromByteArray(byte[] bytes) {
        Roi roi = null;
        try {
            RoiDecoder decoder = new RoiDecoder(bytes, null);
            roi = decoder.getRoi();
        } catch (IOException e) {
            return null;
        }
        return roi;
    }

}
