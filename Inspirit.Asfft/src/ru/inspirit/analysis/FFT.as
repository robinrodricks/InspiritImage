package ru.inspirit.analysis 
{
	import flash.utils.Endian;
	import com.joa_ebert.apparat.memory.Memory;
	import com.joa_ebert.apparat.memory.MemoryMath;
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
	public final class FFT 
	{
		public static const REAL_RAW_DATA:int = 0;
		public static const IMAGINARY_RAW_DATA:int = 1;
		public static const REAL_FFT_DATA:int = 2;
		public static const IMAGINARY_FFT_DATA:int = 3;
		
		public static const TEMP_DATA:int = 7;
		//public static const SPECTRUM_DATA:int = 4;

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

		protected var length:int;
		protected var numChannels:int;
		protected var specLength:int;

		public function FFT():void
		{
			FFT_LIB = (new CLibInit()).init();
			
			var ns:Namespace = new Namespace( "cmodule.fft" );
			alchemyRAM = (ns::gstate).ds;

			getMemoryPointers();
		}
		
		/**
		 * Init FFT instance of specified length
		 * also you may choose to work with STEREO (2 channels) or MONO (1 channel)
		 */
		public function init(length:int = 2048, numChannels:int = 1):int
		{
			this.numChannels = numChannels;
			this.length = MemoryMath.nextPow2(length);
			
			FFT_LIB.initSignalBuffers(this.length, this.numChannels);
			
			specLength = (this.length >> 1) + 1;
			
			return this.length;
		}
		
		/**
		 * Perform forward FFT
		 */
		public function forwardFFT():void
		{
			FFT_LIB.analyzeSignal(1, 0, 0, 0);
		}
		
		/**
		 * Perform inverse FFT
		 */
		public function inverseFFT():void
		{
			FFT_LIB.analyzeSignal(0, 0, 1, 0);
		}
		
		/**
		 * Calculate Spectrum (amgnitude) of Real and Imaginary FFT parts
		 * see getSpectrumData method for more info
		 */
		public function calculateSpectrum(normalizeSpectrum:Boolean = false):void
		{
			FFT_LIB.analyzeSignal(0, 1, 0, normalizeSpectrum ? 1 : 0);
		}
		
		/**
		 * Result Spectrum length
		 */
		public function get spectrumLength():int
		{
			return specLength * numChannels;
		}

		/**
		 * the lib assumes that the input data is formated [left, right, left, right, ...]
		 * the way it is presented in Sound raw data
		 */
		public function setStereoRAWDataByteArray(data:ByteArray):void
		{ 
			alchemyRAM.position = Memory.readInt(shiftedDataPtr);
			alchemyRAM.writeBytes(data);
			
			FFT_LIB.splitChannels();
		}
		
		/**
		 * Use it if you are working with Sound raw data [left, right, left, right, ...]
		 * It will format data before return
		 * May be used after calling inverseFFT to get modified sound signal
		 * Please note data that comes from Alchemy is always LITTLE_ENDIAN
		 */
		public function getStereoRAWDataByteArray(endian:String = Endian.LITTLE_ENDIAN):ByteArray
		{ 
			FFT_LIB.mergeChannels();
			
			var pos:int = Memory.readInt(shiftedDataPtr);
			
			var ba:ByteArray = new ByteArray();
			ba.endian = endian;
			
			ba.writeBytes(alchemyRAM, pos, (length << 1) << 2);
			ba.position = 0;
			
			return ba;
		}
		
		public function setLeftChannelDataByteArray(data:ByteArray, dataType:int = 0):void
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]);
			
			alchemyRAM.position = pos;
			alchemyRAM.writeBytes(data);
		}
		
		public function setRightChannelDataByteArray(data:ByteArray, dataType:int = 0):void
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]) + (length << 2);
			
			alchemyRAM.position = pos;
			alchemyRAM.writeBytes(data);
		}
		
		public function setLeftChannelDataVector(data:Vector.<Number>, dataType:int = 0):void
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]);
			
			var len:int = length;
			var ind1:int = -1;
			var ind2:int = 0;
			
			for(; ind2 < len; ++ind2)
			{
				Memory.writeFloat(data[++ind1], pos + (ind2 << 2));
			}
		}
		
		public function setRightChannelDataVector(data:Vector.<Number>, dataType:int = 0):void
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]) + (length << 2);
			
			var len:int = length;
			var ind1:int = -1;
			var ind2:int = 0;
			
			for(; ind2 < len; ++ind2)
			{
				Memory.writeFloat(data[++ind1], pos + (ind2 << 2));
			}
		}
		
		public function getLeftChannelDataByteArray(dataType:int = 0, endian:String = Endian.LITTLE_ENDIAN):ByteArray
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]);
			
			var ba:ByteArray = new ByteArray();
			ba.endian = endian;
			
			ba.writeBytes(alchemyRAM, pos, length << 2);
			ba.position = 0;
			
			return ba;
		}
		
		public function getRightChannelDataByteArray(dataType:int = 0, endian:String = Endian.LITTLE_ENDIAN):ByteArray
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]) + (length << 2);
			
			var ba:ByteArray = new ByteArray();
			ba.endian = endian;
			
			ba.writeBytes(alchemyRAM, pos, length << 2);
			ba.position = 0;
			
			return ba;
		}
		
		public function getLeftChannelDataVector(dataType:int = 0):Vector.<Number>
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]);
			
			var len:int = length;
			var ind1:int = -1;
			var ind2:int = 0;
			
			var data:Vector.<Number> = new Vector.<Number>(len, true);
			
			for(; ind2 < len; ++ind2)
			{
				data[++ind1] = Memory.readFloat(pos + (ind2 << 2));
			}
			
			return data;
		}
		
		public function getRightChannelDataVector(dataType:int = 0):Vector.<Number>
		{
			var pos:int = Memory.readInt(DATA_POINTERS[dataType]) + (length << 2);
			
			var len:int = length;
			var ind1:int = -1;
			var ind2:int = 0;
			
			var data:Vector.<Number> = new Vector.<Number>(len, true);
			
			for(; ind2 < len; ++ind2)
			{
				data[++ind1] = Memory.readFloat(pos + (ind2 << 2));
			}
			
			return data;
		}
		
		/**
		 * Note: spectum size is (SignalLength / 2 + 1) wide
		 * if you work with 2 channels it will double
		 * half for the left channel and half for the right
		 */
		public function getSpectrumData():ByteArray
		{
			var ba:ByteArray = new ByteArray();
			ba.endian = Endian.LITTLE_ENDIAN;
			
			var pos:int = Memory.readInt(amplFFTDataPtr);
			
			ba.writeBytes(alchemyRAM, pos, (specLength * numChannels) << 2);
			ba.position = 0;
			
			return ba;
		}
		
		/**
		 * Get access to core Alchemy Lib object
		 * needed for FFTSpectrumAnalyzer class
		 */
		public function get core():Object
		{
			return FFT_LIB;
		}
		
		/**
		 * get access to Alchemy Lib Memory ByteArray
		 */
		public function get coreMemory():ByteArray
		{
			return alchemyRAM;
		}

		/**
		 * get memory offset for read/write specified data to/from Alchemy
		 */
		public function getDataMemoryOffset(dataType:int):int
		{
			return Memory.readInt(DATA_POINTERS[dataType]);
		}
		
		/**
		 * return number of channels in use
		 */
		public function get numberOfChannels():int
		{
			return numChannels;
		}
		
		/**
		 * return provided at init length
		 */
		public function get signalLength():int
		{
			return length;
		}

		/**
		 * Clear all allocated buffers
		 * Should be caled before reiniting instance with new settings
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
