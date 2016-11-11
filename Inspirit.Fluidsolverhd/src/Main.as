package 
{
	import ru.inspirit.utils.FluidSolverHD;
	import ru.inspirit.utils.ColorUtils;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.BlendMode;
	import flash.display.Sprite;
	import flash.display.StageScaleMode;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import flash.filters.BlurFilter;
	import flash.geom.ColorTransform;
	import flash.geom.Matrix;
	import flash.geom.Point;
	import flash.ui.ContextMenu;
	import flash.ui.ContextMenuItem;
	import flash.utils.getTimer;

	/**
     * Fluid Solver test application [Alchemy version]
     * @author Eugene Zatepyakin
     * @link http://blog.inspirit.ru/?p=339
     * @link http://code.google.com/p/in-spirit/source/browse/#svn/trunk/projects/FluidSolverHD
     */
	 
	[SWF(width='960',height='580',frameRate='60',backgroundColor='0x000000')]
	
	public class Main extends Sprite 
	{
		
		protected const origin:Point = new Point();
		protected const identity:Matrix = new Matrix();
		protected const blur:BlurFilter = new BlurFilter( 2, 2, 2 );
		protected const fade2alpha:ColorTransform = new ColorTransform( 1, 1, 1, .9);
		
		public var sw:uint = 960;
		public var sh:uint = 540;
		
		public var FLUID_WIDTH:uint = 100;
		
		public var isw:Number = 1 / sw;
		public var ish:Number = 1 / sh;
		
		public var aspectRatio:Number = sw * ish;
		public var aspectRatio2:Number = aspectRatio * aspectRatio;
		
		protected const fSolver:FluidSolverHD = new FluidSolverHD( FLUID_WIDTH, int( FLUID_WIDTH * sh / sw ), sw, sh );
		protected const particlesImage:BitmapData = fSolver.particlesImage;
		protected const particlesBitmap:Bitmap = new Bitmap(new BitmapData(sw, sw, true, 0), 'never', false);
		
		protected var sparkleScale:int = 2; 
		private var sparkle:BitmapData = new BitmapData(sw/sparkleScale, sh/sparkleScale, true, 0x0);
		private var sparkleDrawMatrix:Matrix = new Matrix(1/sparkleScale, 0, 0, 1/sparkleScale, 0, 0);

		protected var prevMouse:Point = new Point();
		protected var frameCount:uint = 0;
		protected var controls:Controls;
		
		protected var drawMode:int = 0;
		protected var mouseDrag:Boolean = false;
		
		public function Main()
		{
			if (stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		protected function init(e:Event = null):void 
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			initStage();
			
			fSolver.width = sw;
			fSolver.height = sh;
			fSolver.y = particlesBitmap.y = 40;			
			
			controls = new Controls(fSolver);
			
			var b:Bitmap = new Bitmap(sparkle, 'never', true);
            b.y = 40;
            b.blendMode = BlendMode.ADD;
            b.scaleX = b.scaleY = sparkleScale;
            
			//particlesBitmap.blendMode = BlendMode.OVERLAY;
			
			addChild(fSolver);
			addChild(particlesBitmap);
			addChild(b);
			addChild(controls);
			
			addEventListener(Event.ENTER_FRAME, render);
			stage.addEventListener(MouseEvent.MOUSE_MOVE, onMouseMove);
			stage.addEventListener(MouseEvent.MOUSE_DOWN, nextMode);
			stage.addEventListener(MouseEvent.MOUSE_UP, onMouseUp);
		}
		
		protected function render(e:Event):void 
		{			
			var t:int = getTimer();
			
			fSolver.update();
			
			if(drawMode > 0)
			{
				particlesBitmap.bitmapData.lock();
				if(drawMode == 1 || drawMode == 3) {
					particlesBitmap.bitmapData.copyPixels(particlesImage, particlesImage.rect, origin);
				} else if(drawMode == 2 || drawMode == 4){
					particlesBitmap.bitmapData.colorTransform(particlesBitmap.bitmapData.rect, fade2alpha);
					particlesBitmap.bitmapData.copyPixels(particlesImage, particlesImage.rect, origin, particlesImage, origin, true);
				}
				particlesBitmap.bitmapData.unlock(particlesBitmap.bitmapData.rect);
				
				sparkle.lock();
	            sparkle.fillRect(sparkle.rect, 0x0);
	            sparkle.draw(particlesBitmap.bitmapData, sparkleDrawMatrix);
	            sparkle.unlock(sparkle.rect);
			}
			
			controls.alg_t += getTimer() - t;
			controls.alg_n ++;
			
			frameCount = ++frameCount % 0xFFFFFFFF;
		}
		
		protected function nextMode(e:MouseEvent = null):void
		{
			if (mouseY <= 40) return;
			
			mouseDrag = true;
			
			drawMode = ++drawMode % 5;
			fSolver.drawMode = drawMode;
			particlesBitmap.bitmapData.fillRect(particlesBitmap.bitmapData.rect, 0);
			sparkle.fillRect(sparkle.rect, 0);
			fSolver.bitmapData.fillRect(fSolver.bitmapData.rect, 0);
		}
		protected function onMouseUp(e:MouseEvent = null):void
		{
			mouseDrag = false;
		}

		protected function onMouseMove(e:MouseEvent):void 
		{
			if (mouseY <= 40) return;
			
			var NormX:Number = mouseX * isw;
			var NormY:Number = (mouseY - 40) * ish;
			var VelX:Number = (mouseX - prevMouse.x) * isw;
			var VelY:Number = (mouseY - 40 - prevMouse.y) * ish;

			addForce(NormX, NormY, VelX, VelY);
			
			prevMouse.x = mouseX;
			prevMouse.y = mouseY-40;
		}
		
		protected function addForce(x:Number, y:Number, dx:Number, dy:Number):void
		{
			var speed:Number = dx * dx  + dy * dy * aspectRatio2;    // balance the x and y components of speed with the screen aspect ratio
			
			if(speed > 0) {
				if (x < 0) x = 0;
				else if (x > 1) x = 1;
				if (y < 0) y = 0;
				else if (y > 1) y = 1;
				
				const hue:Number = ((x + y) * 180 + frameCount) % 360;
				
				fSolver.addForce(x, y, dx, dy, ColorUtils.HSB2GRB(hue), true, true);
			}
		}
		
		private function initStage():void
		{
			stage.scaleMode = StageScaleMode.NO_SCALE;
			//stage.align = StageAlign.TOP_LEFT;
			
			var myContextMenu:ContextMenu = new ContextMenu();
			myContextMenu.hideBuiltInItems();
			
			
			var copyr:ContextMenuItem = new ContextMenuItem("© inspirit.ru", true, false);
			myContextMenu.customItems.push(copyr);
			
			contextMenu = myContextMenu;
		}
	}
}
