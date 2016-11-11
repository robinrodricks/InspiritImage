package demo
{
	import apparat.memory.Memory;

	import ru.inspirit.image.filter.GaussianFilter;

	import com.bit101.components.HUISlider;
	import com.bit101.components.PushButton;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.filters.ColorMatrixFilter;
	import flash.filters.ConvolutionFilter;
	import flash.geom.Point;
	import flash.text.TextField;
	import flash.ui.ContextMenu;
	import flash.ui.ContextMenuItem;
	import flash.utils.ByteArray;
	import flash.utils.Endian;
	import flash.utils.getTimer;

	/**
	 * @author Eugene Zatepyakin
	 */
	[SWF(width='800',height='700',frameRate='25',backgroundColor='0xFFFFFF')]
	public final class GaussTest extends Sprite
	{
		[Embed(source = '../../assets/graffiti_800.jpg')] private static const img_ass:Class;
		
		public const ORIGIN:Point = new Point();
		public static const GAUSSIAN_3x3:ConvolutionFilter = new ConvolutionFilter(
																				3,3,
																				[ 	1,2,1,
																					2,4,2,
																					1,2,1	], 
																				16);
		public static const GAUSSIAN_5x5:ConvolutionFilter = new ConvolutionFilter(
																					5, 5,
																					[	2, 4, 5,4,2,
																						4, 9,12,9,4,
																						5,12,15,12,5,
																						4, 9,12,9,4,
																						2, 4, 5,4,2	],
																						159);
		public static const GRAYSCALE_MATRIX:ColorMatrixFilter = new ColorMatrixFilter([
                        																.2989, .587, .114, 0, 0,
            																			.2989, .587, .114, 0, 0,
            																			.2989, .587, .114, 0, 0,
            																			0, 0, 0, 0, 0
																						]);
																				
		protected var _txt:TextField;

		public const ram:ByteArray = new ByteArray();
		
		public var bmp:BitmapData = Bitmap(new img_ass()).bitmapData;
		public var scr:Bitmap = new Bitmap(bmp);
		
		public var NUM_ITER:int = 1;

		public var frame:int = 0;
		public var frame_prev:int = 0;
		public var frame_ms:int = 0;
		public var aver_time:int = 0;
		public var aver_prev:int = 0;
		public var aver_n:int = 0;
		
		public var imgPtr:int;

		public function GaussTest()
		{
			if(stage)
				init();
			else
				addEventListener( Event.ADDED_TO_STAGE, init );
		}

		protected function init( e:Event = null ):void
		{
			removeEventListener( Event.ADDED_TO_STAGE, init );
			//
			stage.scaleMode = StageScaleMode.NO_SCALE;
			var myContextMenu:ContextMenu = new ContextMenu();
			myContextMenu.hideBuiltInItems();
			var copyr:ContextMenuItem = new ContextMenuItem("© inspirit.ru", true, false);
			myContextMenu.customItems.push(copyr);			
			contextMenu = myContextMenu;

			// debug test field
			_txt = new TextField();
			_txt.autoSize = 'left';
			_txt.width = 300;
			_txt.x = 800 - 140;
			_txt.y = 645;
			addChild( _txt );
			
			//bmp.noise(10, 0, 255, 0, true);
			bmp.applyFilter(bmp, bmp.rect, ORIGIN, GRAYSCALE_MATRIX);
			addChild(scr);

			ram.endian = Endian.LITTLE_ENDIAN;
			ram.length = 1024 + (bmp.width * bmp.height) + ((bmp.width * bmp.height) << 2);
			ram.position = 0;
			Memory.select( ram );
			
			// feed the ram with data
			var data:Vector.<uint> = bmp.getVector( bmp.rect );
			var n:int = data.length;
			var ptr:int = 1024;
			for(var i:int = 0; i < n; ++i, ++ptr)
			{
				Memory.writeByte(data[i], ptr);
			}
			
			imgPtr = 1024;
			
			var pb:PushButton = new PushButton(this, 10, 670, 'CONV FILTER 3x3', testGaussConvFilter3x3);
			new PushButton(this, 120, 670, 'CONV FILTER 5x5', testGaussConvFilter5x5);
			pb = new PushButton(this, 230, 670, 'FMEM GAUSS 3x3', testGaussFastMem3x3);
			//pb.width = 130;
			pb = new PushButton(this, 340, 670, 'FMEM GAUSS 5x5', testGauss5x5);
			pb = new PushButton(this, 450, 670, 'FMEM GAUSS 7x7', testGauss7x7);
			pb = new PushButton(this, 230, 645, 'RESET IMG', resetImage);
			
			var sl:HUISlider = new HUISlider(this, 10, 644, 'NUMBER OF ITER', onIterChange);
			sl.width = 240;
			sl.setSliderParams(1, 20, 1);
			sl.tick = 1;
			
			//stage.addEventListener(MouseEvent.CLICK, runTests);
		}
		
		public function runTests(e:Event = null):void
		{
			_txt.text = '';
			//setTimeout(testGaussConvFilter, 1000);
			//setTimeout(testGaussFastMem, 2000);
			//setTimeout(testGaussFastMemAsm, 3000);
			//testGaussConvFilter();
			//testGaussFastMemAsm();
			//testGauss5x5();
			//testGauss7x7();
			//testGaussFastMem();
		}
		public function resetImage(e:Event = null):void
		{
			bmp = Bitmap(new img_ass()).bitmapData;
			bmp.applyFilter(bmp, bmp.rect, ORIGIN, GRAYSCALE_MATRIX);
			scr.bitmapData = bmp;
			// feed the ram with data
			var data:Vector.<uint> = bmp.getVector( bmp.rect );
			var n:int = data.length;
			var ptr:int = imgPtr;
			for(var i:int = 0; i < n; ++i, ++ptr)
			{
				Memory.writeByte(data[i], ptr);
			}
		}
		
		public function onIterChange(e:Event):void
		{
			NUM_ITER = HUISlider(e.currentTarget).value;
		}
		
		public function testGaussConvFilter3x3(e:Event = null):void
		{
			_txt.text = '';
			var t:int = getTimer();
			for(var i:int = 0; i < NUM_ITER; ++i)
			{
				bmp.applyFilter( bmp, bmp.rect, ORIGIN, GAUSSIAN_3x3 );
			}
			t = (getTimer() - t);
			_txt.appendText('Conv Filter 3x3 \ntotal: ' + (t) + 'ms\n');
			_txt.appendText('average: ' + int(t/NUM_ITER+0.5) + 'ms\n');
		}
		public function testGaussConvFilter5x5(e:Event = null):void
		{
			_txt.text = '';
			var t:int = getTimer();
			for(var i:int = 0; i < NUM_ITER; ++i)
			{
				bmp.applyFilter( bmp, bmp.rect, ORIGIN, GAUSSIAN_5x5 );
			}
			t = (getTimer() - t);
			_txt.appendText('Conv Filter 5x5 \ntotal: ' + (t) + 'ms\n');
			_txt.appendText('average: ' + int(t/NUM_ITER+0.5) + 'ms\n');
		}
		
		public function testGaussFastMem3x3(e:Event = null):void
		{
			var w:int = bmp.width;
			var h:int = bmp.height;
			var ptr:int = 1024;
			var ptri:int = ptr + (w * h); // intermediate image
			_txt.text = '';
			var t:int = getTimer();
			for(var i:int = 0; i < NUM_ITER; ++i)
			{
				GaussianFilter.gaussSmooth3x3Standard(ptr, ptr, w, h, ptri);
			}
			t = (getTimer() - t);
			_txt.appendText('Fast Mem GAUSS 3x3 \ntotal: ' + (t) + 'ms\n');
			_txt.appendText('average: ' + int(t/NUM_ITER+0.5) + 'ms\n');
			renderMem(bmp);
		}
		
		public function testGauss5x5(e:Event = null):void
		{
			var w:int = bmp.width;
			var h:int = bmp.height;
			var ptr:int = 1024;
			var ptri:int = ptr + (w * h); // intermediate image
			_txt.text = '';
			var t:int = getTimer();
			for(var i:int = 0; i < NUM_ITER; ++i)
			{
				GaussianFilter.gaussSmooth5x5Standard(ptr, ptr, w, h, ptri);
			}
			t = (getTimer() - t);
			_txt.appendText('Fast Mem GAUSS 5x5 \ntotal: ' + (t) + 'ms\n');
			_txt.appendText('average: ' + int(t/NUM_ITER+0.5) + 'ms\n');
			renderMem(bmp);
		}
		public function testGauss7x7(e:Event = null):void
		{
			var w:int = bmp.width;
			var h:int = bmp.height;
			var ptr:int = 1024;
			var ptri:int = ptr + (w * h); // intermediate image
			_txt.text = '';
			var t:int = getTimer();
			for(var i:int = 0; i < NUM_ITER; ++i)
			{
				GaussianFilter.gaussSmooth7x7Standard(ptr, ptr, w, h, ptri);
			}
			t = (getTimer() - t);
			_txt.appendText('Fast Mem GAUSS 7x7 \ntotal: ' + (t) + 'ms\n');
			_txt.appendText('average: ' + int(t/NUM_ITER+0.5) + 'ms\n');
			renderMem(bmp);
		}
		
		public function renderMem(bmp:BitmapData):void
		{
			var w:int = bmp.width;
			var h:int = bmp.height;
			var area:int = w * h;
			var vec:Vector.<uint> = new Vector.<uint>( area );
			var _p:int = 1024;
			for(var i:int = 0; i < area; ++i, ++_p)
			{
				var c:uint = Memory.readUnsignedByte(_p);
				vec[i] = c << 16 | c << 8 | c;
			}
			bmp.lock();
			bmp.setVector(bmp.rect, vec);
			bmp.unlock();
		}
	}
}
