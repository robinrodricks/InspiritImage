package ru.inspirit.analysis
{
	import flash.geom.Rectangle;
	import flash.utils.Endian;
	import com.joa_ebert.apparat.memory.MemoryMath;
	import com.joa_ebert.apparat.memory.Memory;
	import flash.display.BitmapData;
	import cmodule.fft.CLibInit;

	import flash.utils.ByteArray;

	/**
	*  released under MIT License (X11)
	*  http://www.opensource.org/licenses/mit-license.php
	*
	*  Eugene Zatepyakin
	*  http://blog.inspirit.ru
	*  http://code.google.com/p/in-spirit/wiki/ASFFT
	*/
	public final class FFT2D
	{
		public static const REAL_RAW_DATA:int = 0;
		public static const IMAGINARY_RAW_DATA:int = 1;
		public static const REAL_FFT_DATA:int = 2;
		public static const IMAGINARY_FFT_DATA:int = 3;
		public static const AMPLITUDE_DATA:int = 4;
		public static const PHASE_DATA:int = 5;

		protected const DATA_POINTERS:Vector.<int> = new Vector.<int>(8, true);
		
		protected var FFT_LIB:Object;
		protected var alchemyRAM:ByteArray;

		protected var realDataPtr:int;
		protected var imagDataPtr:int;
		protected var realFFTDataPtr:int;
		protected var imagFFTDataPtr:int;
		protected var amplFFTDataPtr:int;
		protected var phaseFFTDataPtr:int;
		protected var drawDataPtr:int;
		protected var shiftedDataPtr:int;

		protected var imageW:int;
		protected var imageH:int;
		protected var imageW2:int;
		protected var imageH2:int;
		protected var area2:int;
		protected var numChannels:int;

		public function FFT2D():void
		{
			FFT_LIB = (new CLibInit()).init();
			
			var ns:Namespace = new Namespace( "cmodule.fft" );
			alchemyRAM = (ns::gstate).ds;

			getMemoryPointers();
		}
		
		/**
		 * Returned Rectangle represents actual data dimensions
		 * inside Lib. width and height rounded to min/max power of 2
		 */
		public function init(width:int, height:int, numChannels:int = 3, highPrecision:Boolean = true):Rectangle
		{
			var w:int = imageW = width;
			var h:int = imageH = height;
			
			var w2:int;
			var h2:int;
			
			if(highPrecision)
			{
				w2 = imageW2 = Math.max(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			} else {
				w2 = imageW2 = Math.min(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			}
			
			h2 = imageH2 = w2;
			
			this.numChannels = numChannels;
			area2 = imageW2 * imageH2 * numChannels;

			FFT_LIB.allocateBuffers(w, h, w2, h2, numChannels);
			
			return new Rectangle(0, 0, w2, h2);
		}
		
		/**
		 * All 3 Channels used as inut
		 * data buffer size is extended to next power of 2 (if needed)
		 */
		public function initFromRGBBitmap(bmp:BitmapData, highPrecision:Boolean = true):void
		{
			var w:int = imageW = bmp.width;
			var h:int = imageH = bmp.height;

			var w2:int;
			var h2:int;
			
			if(highPrecision)
			{
				w2 = imageW2 = Math.max(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			} else {
				w2 = imageW2 = Math.min(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			}
			
			h2 = imageH2 = w2;

			area2 = (imageW2 * imageH2 * 3);
			numChannels = 3;

			FFT_LIB.allocateBuffers(w, h, w2, h2, numChannels);

			var realPos:int = Memory.readInt(realDataPtr);

			var c:uint;
			var ind:int;

			for(var i:int = 0; i < h; ++i)
			{
				ind = (i * w2) * 3;
				for(var j:int = 0; j < w; ++j, ++ind)
				{
					c = bmp.getPixel(j, i);

					Memory.writeFloat(c >> 16 & 0xFF, realPos + (ind << 2));
					Memory.writeFloat(c >> 8 & 0xFF, realPos + ((++ind) << 2));
					Memory.writeFloat(c & 0xFF, realPos + ((++ind) << 2));
				}
			}
		}

		/**
		 * Only BLUE Channel data is used as input
		 * data buffer size is extended to next power of 2 (if needed)
		 */
		public function initFromGrayBitmap(bmp:BitmapData, highPrecision:Boolean = true):void
		{
			var w:int = imageW = bmp.width;
			var h:int = imageH = bmp.height;

			var w2:int;
			var h2:int;
			
			if(highPrecision)
			{
				w2 = imageW2 = Math.max(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			} else {
				w2 = imageW2 = Math.min(MemoryMath.nextPow2(w), MemoryMath.nextPow2(h));
			}
			
			h2 = imageH2 = w2;

			area2 = imageW2 * imageH2;
			numChannels = 1;

			FFT_LIB.allocateBuffers(w, h, w2, h2, numChannels);

			var realPos:int = Memory.readInt(realDataPtr);

			var ind:int;

			for(var i:int = 0; i < h; ++i)
			{
				ind = (i * w2);
				for(var j:int = 0; j < w; ++j, ++ind)
				{

					Memory.writeFloat(bmp.getPixel(j, i) & 0xFF, realPos + (ind << 2));
				}
			}
		}
		
		/**
		 * Perform forward FFT
		 */
		public function forwardFFT():void
		{
			FFT_LIB.analyzeImage(1, 0, 0, 0);
		}

		/**
		 * Perform inverse FFT
		 */
		public function inverseFFT():void
		{
			FFT_LIB.analyzeImage(0, 0, 0, 1);
		}
	
		/**
		 * Calculate amplitude of Real and Imaginary FFT parts
		 * sqrt(real*real + imaginary*imaginary);  
		 */
		public function calculateAmplitude():void
		{
			FFT_LIB.analyzeImage(0, 1, 0, 0);
		}
		
		/**
		 * Calculate phase of Real and Imaginary FFT parts
		 * atan2(real, imaginary)
		 */
		public function calculatePhase():void
		{
			FFT_LIB.analyzeImage(0, 0, 1, 0);
		}

		/**
		 * Short cut to perform several actions
		 */
		public function analyzeImage(getFFT:Boolean = true, getAmplitude:Boolean = false, getPhase:Boolean = false, getIFFT:Boolean = false):void
		{
			FFT_LIB.analyzeImage(getFFT ? 1 : 0, getAmplitude ? 1 : 0, getPhase ? 1 : 0, getIFFT ? 1 : 0);
		}
		
		/**
		 * Set specified data from BitmapData
		 */
		public function setDataBitmapData(bmp:BitmapData, dataType:int = 0):void
		{
			var w:int= bmp.width;
			var h:int= bmp.height;
			var w2:int = imageW2;
			var h2:int = imageH2;
			var hw:int = w2 >> 1;
			var hh:int = h2 >> 1;
			
			var maxd:Number = Math.sqrt( 2 * hw * hh );
			var d:Number;
			
			var ind:int;
			var i:int, j:int, x:int, y:int;
			var c:uint;
			var pix:int;
			
			var realPos:int = Memory.readInt(DATA_POINTERS[dataType]);
			
			if(numChannels == 3) 
			{
				for(i = 0; i < h; ++i)
				{
					ind = (i * w2) * 3;
					for(j = 0; j < w; ++j, ++ind)
					{
						c = bmp.getPixel(j, i);
						
						// shift corners
						ind = (((i + hh) % h2) * w2 + ((j + hw) % w2)) * 3;
						
						x = j-hw;
			 			y = i-hh;
						
						d = maxd / ( 1 +  Math.sqrt( x*x+y*y ) );
						
						pix = c >> 16 & 0xFF;
						
						Memory.writeFloat(((Math.abs(pix-127)*1.7*d) * ( pix < 127 ? -1 : 1)), realPos + (ind << 2));
						pix = c >> 8 & 0xFF;
						Memory.writeFloat(((Math.abs(pix-127)*1.7*d) * ( pix < 127 ? -1 : 1)), realPos + ((++ind) << 2));
						pix = c & 0xFF;
						Memory.writeFloat(((Math.abs(pix-127)*1.7*d) * ( pix < 127 ? -1 : 1)), realPos + ((++ind) << 2));
					}
				}
			}
			else
			{
				for(i = 0; i < h; ++i)
				{
					ind = i * w2;
					for(j = 0; j < w; ++j, ++ind)
					{
						// shift corners
						ind = ((i + hh) % h2) * w2 + ((j + hw) % w2);
						
						x = j-hw;
			 			y = i-hh;
						
						d = maxd / ( 1 +  Math.sqrt( x*x+y*y ) );
						
						pix = bmp.getPixel(j, i) & 0xFF;
						
						Memory.writeFloat(((Math.abs(pix-127)*1.7*d) * ( pix < 127 ? -1 : 1)), realPos + (ind << 2));
					}
				}
			}
		}
		
		/**
		 * Get/draw specified data to BitmapData
		 */
		public function getDataBitmapData(bmp:BitmapData, dataType:int = 0):void
		{
			FFT_LIB.drawImagePreserveData(dataType, 1);

			alchemyRAM.position = Memory.readInt(drawDataPtr);

			bmp.lock();
			bmp.setPixels(bmp.rect, alchemyRAM);
			bmp.unlock();
		}

		/**
		 * Set specified data from ByteArray
		 * Note that bytearray should be LITTLE_ENDIAN
		 */
		public function setDataByteArray(data:ByteArray, dataType:int = 0, shiftCorners:Boolean = false):void
		{
			data.position = 0;

			if(shiftCorners)
			{
				alchemyRAM.position = Memory.readInt(shiftedDataPtr);
				alchemyRAM.writeBytes(data);

				FFT_LIB.shiftImageData(dataType, 1);
			} else
			{
				alchemyRAM.position = Memory.readInt(DATA_POINTERS[dataType]);
				alchemyRAM.writeBytes(data);
			}

		}

		/**
		 * Get specified data to ByteArray
		 * Note that returned bytearray would be LITTLE_ENDIAN
		 */
		public function getDataByteArray(dataType:int = 0, shiftCorners:Boolean = false):ByteArray
		{
			var realPos:int;
			var data:ByteArray = new ByteArray();
			data.endian = Endian.LITTLE_ENDIAN;

			if(shiftCorners)
			{
				FFT_LIB.shiftImageData(dataType);
				realPos = Memory.readInt(shiftedDataPtr);

				data.writeBytes(alchemyRAM, realPos, area2 << 2);
			}
			else
			{
				realPos = Memory.readInt(DATA_POINTERS[dataType]);

				data.writeBytes(alchemyRAM, realPos, area2 << 2);
			}

			data.position = 0;

			return data;
		}

		/**
		 * Set specified data from Vector.<Number>
		 */
		public function setDataVector(data:Vector.<Number>, dataType:int = 0, shiftCorners:Boolean = false):void
		{
			var realPos:int = Memory.readInt(DATA_POINTERS[dataType]);

			var dind:int = -1;
			var ind:int = 0;

			if(shiftCorners)
			{
				var w2:int = imageW2;
				var h2:int = imageH2;

				var hw:int = w2 >> 1;
				var hh:int = h2 >> 1;
				var i:int, j:int;

				if(numChannels == 3)
				{
					for(i = 0; i < h2; ++i)
					{
						for(j = 0; j < w2; ++j)
						{
							ind = (((i + hh) % h2) * w2 + ((j + hw) % w2)) * 3;

							Memory.writeFloat(data[++dind], realPos + (ind << 2));
							Memory.writeFloat(data[++dind], realPos + ((++ind) << 2));
							Memory.writeFloat(data[++dind], realPos + ((++ind) << 2));
						}
					}
				} else {
					for(i = 0; i < h2; ++i)
					{
						for(j = 0; j < w2; ++j)
						{
							ind = (((i + hh) % h2) * w2 + ((j + hw) % w2));

							Memory.writeFloat(data[++dind], realPos + (ind << 2));
						}
					}
				}
			}
			else
			{
				if(numChannels == 3)
				{
					for(; ind < area2; ++ind)
					{
						Memory.writeFloat(data[++dind], realPos + (ind << 2));
						Memory.writeFloat(data[++dind], realPos + ((++ind) << 2));
						Memory.writeFloat(data[++dind], realPos + ((++ind) << 2));
					}
				} else {
					for(; ind < area2; ++ind)
					{
						Memory.writeFloat(data[++dind], realPos + (ind << 2));
					}
				}
			}
		}

		/**
		 * Get specified data to Vector.<Number>
		 */
		public function getDataVector(dataType:int = 0, shiftCorners:Boolean = false):Vector.<Number>
		{
			var realPos:int = Memory.readInt(DATA_POINTERS[dataType]);

			var data:Vector.<Number> = new Vector.<Number>(area2, true);
			var dind:int = -1;
			var ind:int = 0;

			if(shiftCorners)
			{
				var w2:int = imageW2;
				var h2:int = imageH2;

				var hw:int = w2 >> 1;
				var hh:int = h2 >> 1;
				var i:int, j:int;
				if(numChannels == 3)
				{
					for(i = 0; i < h2; ++i)
					{
						for(j = 0; j < w2; ++j)
						{
							ind = (((i + hh) % h2) * w2 + ((j + hw) % w2)) * 3;

							data[++dind] = Memory.readFloat(realPos + (ind << 2));
							data[++dind] = Memory.readFloat(realPos + ((++ind) << 2));
							data[++dind] = Memory.readFloat(realPos + ((++ind) << 2));
						}
					}
				} else {
					for(i = 0; i < h2; ++i)
					{
						for(j = 0; j < w2; ++j)
						{
							ind = (((i + hh) % h2) * w2 + ((j + hw) % w2));

							data[++dind] = Memory.readFloat(realPos + (ind << 2));
						}
					}
				}
			}
			else
			{
				if(numChannels == 3)
				{
					for(; ind < area2; ++ind)
					{
						data[++dind] = Memory.readFloat(realPos + (ind << 2));
						data[++dind] = Memory.readFloat(realPos + ((++ind) << 2));
						data[++dind] = Memory.readFloat(realPos + ((++ind) << 2));
					}
				} else {
					for(; ind < area2; ++ind)
					{
						data[++dind] = Memory.readFloat(realPos + (ind << 2));
					}
				}
			}

			return data;
		}

		/**
		 * Draw specified data to bitmapData
		 * Note: shiftCorners is used to see radial (from center to borders) data visual representation
		 * add128 is used with FFT Real and Imaginary data types cause they can be negative
		 */
		public function draw(dataType:int, bmp:BitmapData, shiftCorners:Boolean = false, add128:Boolean = false):void
		{
			var w:int = bmp.width;
			var h:int = bmp.height;

			FFT_LIB.drawImageData(dataType, w, h, shiftCorners ? 1 : 0, add128 ? 1 : 0);

			alchemyRAM.position = Memory.readInt(drawDataPtr);

			bmp.lock();
			bmp.setPixels(bmp.rect, alchemyRAM);
			bmp.unlock();

			/*var w2:int = imageW2;
			var h2:int = imageH2;

			var r:int, g:int, b:int;
			var ind:int = 0;

			var hw:int = w2 >> 1;
			var hh:int = h2 >> 1;
			var stride:int = (w2 - w) * 3;

			bmp.lock();

			for(var i:int = 0; i < h; ++i)
			{
				for(var j:int = 0; j < w; ++j)
				{
					if(shiftCorners)
					{
						ind = (((i + hh) % h2) * w2 + ((j + hw) % w2)) * 3;
					}

					r = data[ind++];
					g = data[ind++];
					b = data[ind++];

					if(add128)
					{
						r += 128;
						g += 128;
						b += 128;
					}

					r = r < 0 ? 0 : (r > 255 ? 255 : r);
					g = g < 0 ? 0 : (g > 255 ? 255 : g);
					b = b < 0 ? 0 : (b > 255 ? 255 : b);

					bmp.setPixel(j, i, r<<16 | g<<8 | b);
				}
				ind += stride;
			}

			bmp.unlock();*/
		}
	
		/**
		 * Clear instance memory buffers
		 * should be performed before reiniting or when decided to finish working with instance 
		 */
		public function clear():void
		{
			FFT_LIB.freeBuffers();
		}

		protected function getMemoryPointers():void
		{
			var ptrs:Array = FFT_LIB.getBufferPointers();

			DATA_POINTERS[0] = realDataPtr = ptrs[0];
			DATA_POINTERS[1] = imagDataPtr = ptrs[1];
			DATA_POINTERS[2] = realFFTDataPtr = ptrs[2];
			DATA_POINTERS[3] = imagFFTDataPtr = ptrs[3];
			DATA_POINTERS[4] = amplFFTDataPtr = ptrs[4];
			DATA_POINTERS[5] = phaseFFTDataPtr = ptrs[7];
			DATA_POINTERS[6] = drawDataPtr = ptrs[5];
			DATA_POINTERS[7] = shiftedDataPtr = ptrs[6];
		}
	}
}
