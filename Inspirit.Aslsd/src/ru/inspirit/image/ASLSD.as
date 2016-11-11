package ru.inspirit.image 
{
	import apparat.memory.Memory;

	import cmodule.lsd.CLibInit;

	import flash.display.BitmapData;
	import flash.filters.BlurFilter;
	import flash.filters.ColorMatrixFilter;
	import flash.geom.Matrix;
	import flash.geom.Point;
	import flash.geom.Rectangle;

	/**
	 * @author Eugene Zatepyakin
	 */
	final public class ASLSD 
	{
		protected static const LSD_LIB:Object = (new CLibInit()).init();
		//protected static const ns:Namespace = new Namespace( "cmodule.lsd" );
		//protected static const alchemyRAM:ByteArray = (ns::gstate).ds;
		
		protected const GRAYSCALE_MATRIX:ColorMatrixFilter = new ColorMatrixFilter([
			.2989, .587, .114, 0, 0,
            .2989, .587, .114, 0, 0,
            .2989, .587, .114, 0, 0,
            0, 0, 0, 0, 0
		]);
		protected const BLUR:BlurFilter = new BlurFilter(2, 2, 2);
		protected const ORIGIN:Point = new Point();
		
		protected var image:BitmapData;
		protected var buffer:BitmapData;
		protected var rect:Rectangle;
		protected var width:int;
		protected var height:int;
		protected var area:int;
		
		protected var inPtr:int;
		protected var outPtr:int;
		protected var cntPtr:int;
		
		protected var quantPtr:int;
		protected var ang_thPtr:int;
		protected var epsPtr:int;
		protected var density_thPtr:int;
		protected var n_binsPtr:int;
		protected var max_gradPtr:int;
		
		protected var scale:Number = 0.8;
		protected var iscale:Number = 1 / scale;
		protected var scaleMatrix:Matrix = new Matrix(scale, 0, 0, scale);

		public var segmentsCount:int = 0;
		
		public function ASLSD(downSample:Number = 0.8, quant:Number = 2.0, gradientAngleTolerance:Number = 22.5, 
								detectionThreshold:Number = 0.0, densityOfRegionPoints:Number = 0.7, gradientBins:int = 1024)
		{
			this.downSample = downSample;
			this.quant = quant;
			this.gradientAngleTolerance = gradientAngleTolerance;
			this.detectionThreshold = detectionThreshold;
			this.densityOfRegionPoints = densityOfRegionPoints;
			this.gradientBins = gradientBins;
		}
		
		public function getLineSegments():Vector.<int>
		{
			if(scale < 1.0) 
			{
				buffer.draw(image, scaleMatrix);
				buffer.applyFilter(buffer, rect, ORIGIN, GRAYSCALE_MATRIX);
			} 
			else 
			{
				buffer.applyFilter(image, rect, ORIGIN, GRAYSCALE_MATRIX);
			}
			
			if(blurSize > 0) buffer.applyFilter(buffer, rect, ORIGIN, BLUR);
			
			var data:Vector.<uint> = buffer.getVector(rect);
			var i:int = area;
			while( --i > -1 )
			{
				Memory.writeDouble((data[i] & 0xFF), inPtr + (i << 3));
			}
			
			LSD_LIB.runLSD();
			
			i = segmentsCount = Memory.readInt(cntPtr);
			var outp:int = Memory.readInt(outPtr);
			
			var res:Vector.<int> = new Vector.<int>(i << 2, true);
			var j:int = -1;
			while( --i > -1 )
			{
				res[++j] = Memory.readDouble(outp + ((i*5) << 3)) * iscale;
				res[++j] = Memory.readDouble(outp + ((i*5+1) << 3)) * iscale;
				res[++j] = Memory.readDouble(outp + ((i*5+2) << 3)) * iscale;
				res[++j] = Memory.readDouble(outp + ((i*5+3) << 3)) * iscale;
			}
			
			return res;
		}
		
		public function dispose():void
		{
			LSD_LIB.dispose();
			buffer.dispose();
		}

		public function set source(bmp:BitmapData):void
		{
			this.image = bmp;
			
			if(bmp.width * scale != width || bmp.height * scale != height)
			{
				downSample = scale;
			}
		}
		
		public function set downSample(scale:Number):void
		{
			this.scale = scale;
			this.iscale = 1 / scale;
			
			scaleMatrix.identity();
			scaleMatrix.scale(scale, scale);
			
			if(this.image)
			{
				width = this.image.width * scale;
				height = this.image.height * scale;
				area = width * height;
				
				buffer = new BitmapData(width, height, false, 0x00);
				buffer.lock();
				
				rect = buffer.rect;
				
				var pts:Array = LSD_LIB.setupGlobalBuffers(width, height);
				
				inPtr = pts[0];
				outPtr = pts[1];
				cntPtr = pts[2];
				
				quantPtr = pts[3];
				ang_thPtr = pts[4];
				epsPtr = pts[5];
				density_thPtr = pts[6];
				n_binsPtr = pts[7];
				max_gradPtr = pts[8];
			}
		}
		
		public function set blurSize(val:int):void
		{
			BLUR.blurX = BLUR.blurY = val;
		}
		public function get blurSize():int
		{
			return BLUR.blurX;
		}

		/**
		 *  Bound to the quantization error on the gradient norm.
		 *  Example: if gray level is quantized to integer steps,
		 *  the gradient (computed by finite differences) error
		 *  due to quantization will be bounded by 2.0, as the
		 *  worst case is when the error are 1 and -1, that
		 *  gives an error of 2.0.
		 *  Suggested value: 2.0
		 */
		public function set quant(val:Number):void
		{
			Memory.writeDouble(val, quantPtr);
		}
		public function get quant():Number
		{
			return Memory.readDouble(quantPtr);
		}
		
		/**
		 * Gradient angle tolerance in the region growing
		 * algorithm, in degrees.
		 * Suggested value: 22.5
		 */
		public function set gradientAngleTolerance(val:Number):void
		{
			Memory.writeDouble(val, ang_thPtr);
		}
		public function get gradientAngleTolerance():Number
		{
			return Memory.readDouble(ang_thPtr);
		}
		
		/**
		 * Detection threshold, -log10(NFA).
		 * The bigger, the more strict the detector is,
		 * and will result in less detections.
		 * (Note that the 'minus sign' makes that this
		 * behavior is opposite to the one of NFA.)
		 * The value -log10(NFA) is equivalent but more
		 * intuitive than NFA:
		 * 		-1.0 corresponds to 10 mean false alarms
		 * 		0.0 corresponds to 1 mean false alarm
		 * 		1.0 corresponds to 0.1 mean false alarms
		 * 		2.0 corresponds to 0.01 mean false alarms
		 * Suggested value: 0.0
		 */
		public function set detectionThreshold(val:Number):void
		{
			Memory.writeDouble(val, epsPtr);
		}
		public function get detectionThreshold():Number
		{
			return Memory.readDouble(epsPtr);
		}
		
		/**
		 * Minimal proportion of region points in a rectangle.
		 * Suggested value: 0.7
		 */
		public function set densityOfRegionPoints(val:Number):void
		{
			Memory.writeDouble(val, density_thPtr);
		}
		public function get densityOfRegionPoints():Number
		{
			return Memory.readDouble(density_thPtr);
		}
		
		/**
		 * Number of bins used in the pseudo-ordering of gradient
		 * modulus.
		 * Suggested value: 1024
		 */
		public function set gradientBins(val:int):void
		{
			Memory.writeInt(val, n_binsPtr);
		}
		public function get gradientBins():int
		{
			return Memory.readInt(n_binsPtr);
		}
	}
}
