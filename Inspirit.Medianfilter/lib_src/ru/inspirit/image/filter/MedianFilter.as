package ru.inspirit.image.filter
{
	import cmodule.median.CLibInit;

	import com.joa_ebert.apparat.memory.Memory;

	import flash.display.BitmapData;
	import flash.display.Shader;
	import flash.display.ShaderJob;
	import flash.events.Event;
	import flash.events.EventDispatcher;
	import flash.events.ProgressEvent;
	import flash.events.TimerEvent;
	import flash.utils.ByteArray;
	import flash.utils.Timer;

	/**
	 * Alchemy based version of Median Filter
	 * 
	 * also included Pixel Bender filters for 3x3 & 5x5 median processing
	 * optimized by Mario Klingemann (http://www.quasimondo.com)
	 *
	 * @author Eugene Zatepyakin
	 * @link http://blog.inspirit.ru
	 */
	final public class MedianFilter extends EventDispatcher
	{
		[Embed('../../../../../pbj/median3x3.pbj', mimeType='application/octet-stream')] private static const median3x3PBJ:Class;
		[Embed('../../../../../pbj/median5x5.pbj', mimeType='application/octet-stream')] private static const median5x5PBJ:Class;

		protected var MEDIAN_LIB:Object;
		protected var alchemyRAM:ByteArray;

		protected var outPointer:int;
		protected var inPointer:int;
		protected var progressPointer:int;

		protected var asyncTimer:Timer;
		protected var buffer:BitmapData;

		protected var median3x3_shader:Shader;		
		protected var median5x5_shader:Shader;
		protected var median_job:ShaderJob;

		public function MedianFilter()
		{
			init();

			median3x3_shader = new Shader(new median3x3PBJ() as ByteArray);
			median5x5_shader = new Shader(new median5x5PBJ() as ByteArray);

			asyncTimer = new Timer(10);
		}

		public function init():void
		{
			MEDIAN_LIB = (new CLibInit()).init();
			var ns:Namespace = new Namespace( "cmodule.median" );
			alchemyRAM = (ns::gstate).ds;

			var points:Array = MEDIAN_LIB.getPointers();

			outPointer = points[0];
			inPointer = points[1];
			progressPointer = points[2];
		}

		public function median3x3(src:BitmapData, dst:BitmapData):void
		{
			median3x3_shader.data.src.width = src.width;
			median3x3_shader.data.src.height = src.height;
			median3x3_shader.data.src.input = src;

			median_job = new ShaderJob(median3x3_shader, dst, dst.width, dst.height);
			median_job.start(true);
		}
		
		public function median5x5(src:BitmapData, dst:BitmapData):void
		{
			median5x5_shader.data.src.width = src.width;
			median5x5_shader.data.src.height = src.height;
			median5x5_shader.data.src.input = src;

			median_job = new ShaderJob(median5x5_shader, dst, dst.width, dst.height);
			median_job.start(true);
		}

		public function medianRGB(src:BitmapData, dst:BitmapData, radius:int):void
		{
			MEDIAN_LIB.setupMedianFilter(src.width, src.height, 3, radius);

			var tmp:int;
			var outp:int = Memory.readInt(outPointer);
			var inp:int = tmp = Memory.readInt(inPointer);

			src.lock();
			var data:Vector.<uint> = src.getVector(src.rect);
			data.fixed = true;

			var len:int = data.length;
			var i:int;

			for(i = 0; i < len; ++i, ++inp)
			{
				Memory.writeInt(data[i], inp++);
				inp++;
			}

			MEDIAN_LIB.runMedianFilter(radius);

			inp = tmp;
			for(i = 0; i < len; ++i, ++outp)
			{
				//data[i] = Memory.readInt(outp++);
				Memory.writeInt(Memory.readInt(outp++), tmp + (i<<2));
				outp++;
			}

			dst.lock();
			alchemyRAM.position = inp;
			dst.setPixels(dst.rect, alchemyRAM);
			//dst.setVector(dst.rect, data);
			dst.unlock(dst.rect);
			src.unlock(src.rect);
		}

		public function medianRGBAsync(src:BitmapData, dst:BitmapData, radius:int):void
		{
			MEDIAN_LIB.setupMedianFilter(src.width, src.height, 3, radius);

			var inp:int = Memory.readInt(inPointer);

			src.lock();
			var data:Vector.<uint> = src.getVector(src.rect);

			var len:int = data.length;
			var i:int;

			for(i = 0; i < len; ++i, ++inp)
			{
				Memory.writeInt(data[i], inp++);
				inp++;
			}

			src.unlock(src.rect);

			buffer = dst;

			MEDIAN_LIB.runMedianFilterAsync(onRGBAsyncComplete, radius);

			asyncTimer.addEventListener(TimerEvent.TIMER, dispatchProgress);
	        asyncTimer.start();
		}

		public function medianGray(src:BitmapData, dst:BitmapData, radius:int):void
		{
			MEDIAN_LIB.setupMedianFilter(src.width, src.height, 1, radius);

			var tmp:int;
			var outp:int = Memory.readInt(outPointer);
			var inp:int = tmp = Memory.readInt(inPointer);

			src.lock();
			var data:Vector.<uint> = src.getVector(src.rect);

			var len:int = data.length;
			var c:uint;
			var i:int;

			for(i = 0; i < len; ++i, ++inp)
			{
				Memory.writeByte(data[i] & 0xFF, inp);
			}

			MEDIAN_LIB.runMedianFilter(radius);

			inp = tmp;
			for(i = 0; i < len; ++i, ++outp)
			{
				c = Memory.readUnsignedByte(outp);
				//data[i] = c << 16 | c << 8 | c;
				Memory.writeInt(c << 16 | c << 8 | c, inp + (i<<2));
			}

			dst.lock();
			alchemyRAM.position = inp;
			dst.setPixels(dst.rect, alchemyRAM);
			//dst.setVector(dst.rect, data);
			dst.unlock(dst.rect);
			src.unlock(src.rect);
		}

		public function medianGrayAsync(src:BitmapData, dst:BitmapData, radius:int):void
		{
			MEDIAN_LIB.setupMedianFilter(src.width, src.height, 1, radius);

			var inp:int = Memory.readInt(inPointer);

			src.lock();
			var data:Vector.<uint> = src.getVector(src.rect);

			var len:int = data.length;
			var i:int;

			for(i = 0; i < len; ++i, ++inp)
			{
				Memory.writeByte(data[i] & 0xFF, inp);
			}

			buffer = dst;

			MEDIAN_LIB.runMedianFilterAsync(onGrayAsyncComplete, radius);

			asyncTimer.addEventListener(TimerEvent.TIMER, dispatchProgress);
	        asyncTimer.start();

			src.unlock(src.rect);
		}

		public function dispose():void
		{
			MEDIAN_LIB.dispose();
			MEDIAN_LIB = null;
		}

		protected function onRGBAsyncComplete(...arg):void
		{
			var outp:int = Memory.readInt(outPointer);
			var inp:int = Memory.readInt(inPointer);
			var i:int;
			var len:int = buffer.width * buffer.height;

			for(i = 0; i < len; ++i, ++outp)
			{
				Memory.writeInt(Memory.readInt(outp++), inp + (i<<2));
				outp++;
			}

			buffer.lock();
			alchemyRAM.position = inp;
			buffer.setPixels(buffer.rect, alchemyRAM);
			buffer.unlock(buffer.rect);

			asyncTimer.stop();
			asyncTimer.removeEventListener(TimerEvent.TIMER, dispatchProgress);

			dispatchEvent(new ProgressEvent(ProgressEvent.PROGRESS, false, false, buffer.width, buffer.width));
			dispatchEvent(new Event(Event.COMPLETE));
		}

		protected function onGrayAsyncComplete(...arg):void
		{
			var outp:int = Memory.readInt(outPointer);
			var inp:int = Memory.readInt(inPointer);

			//var data:Vector.<uint> = buffer.getVector(buffer.rect);
			var len:int = buffer.width * buffer.height;//data.length;

			var c:uint;
			var i:int;

			for(i = 0; i < len; ++i, ++outp)
			{
				c = Memory.readUnsignedByte(outp);
				//data[i] = c << 16 | c << 8 | c;
				Memory.writeInt(c << 16 | c << 8 | c, inp + (i<<2));
			}

			buffer.lock();
			//buffer.setVector(buffer.rect, data);
			alchemyRAM.position = inp;
			buffer.setPixels(buffer.rect, alchemyRAM);
			buffer.unlock(buffer.rect);

			asyncTimer.stop();
			asyncTimer.removeEventListener(TimerEvent.TIMER, dispatchProgress);

			dispatchEvent(new ProgressEvent(ProgressEvent.PROGRESS, false, false, buffer.width, buffer.width));
			dispatchEvent(new Event(Event.COMPLETE));
		}

		protected function dispatchProgress(e:Event):void
		{
			var pc:int = Memory.readInt(progressPointer);
			dispatchEvent(new ProgressEvent(ProgressEvent.PROGRESS, false, false, pc, buffer.width));
		}
	}
}
