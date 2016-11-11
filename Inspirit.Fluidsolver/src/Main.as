package 
{
	import flash.utils.getTimer;
	import flash.geom.Rectangle;
	import ru.inspirit.utils.ColorUtils;
	import ru.inspirit.utils.FluidSolver;

	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.display.BlendMode;
	import flash.display.Sprite;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import flash.filters.BlurFilter;
	import flash.geom.ColorTransform;
	import flash.geom.Matrix;
	import flash.geom.Point;
	
	import flash.ui.ContextMenu;
	import flash.ui.ContextMenuItem;
	import flash.ui.ContextMenuBuiltInItems;
	import flash.display.StageAlign;
	import flash.display.StageScaleMode;

	/**
	 * Fluid Solver test application
	 * @author Eugene Zatepyakin
	 * @link http://blog.inspirit.ru/?p=248
	 * @link http://code.google.com/p/in-spirit/source/browse/#svn/trunk/projects/FluidSolver
	 */
	public class Main extends Sprite 
	{
		private const origin:Point = new Point();
		private const identity:Matrix = new Matrix();
		private const blur:BlurFilter = new BlurFilter( 2, 2, 2 );
		private const fade2black:ColorTransform = new ColorTransform( 0.9, 0.9, 0.9 );
		private const fade2alpha:ColorTransform = new ColorTransform( 1, 1, 1, .7);
		
		public static const sw:uint = 800;
		public static const sh:uint = 600;
		
		private static const DRAW_SCALE:Number = 0.5;
		
		public static const FLUID_WIDTH:uint = 50;
		
		public static const isw:Number = 1 / sw;
		public static const ish:Number = 1 / sh;
		
		public static const aspectRatio:Number = sw * ish;
		public static const aspectRatio2:Number = aspectRatio * aspectRatio;
		
		public static const fSolver:FluidSolver = new FluidSolver( FLUID_WIDTH, int( FLUID_WIDTH * sh / sw ) );
		
		public static var drawFluid:Boolean = true;
		public static var drawParticles:Boolean = true;
		public static var drawLines:Boolean = false;
		public static var drawMode:uint = 0;
		
		public static var screen:BitmapData = new BitmapData(sw, sh, true, 0);
		public static var fluid:BitmapData = new BitmapData(fSolver.width - 2, fSolver.height - 2, false, 0);
		public static var fluidBuffer:Vector.<uint> = new Vector.<uint>(fluid.width * fluid.height, true);
		
		private var fade:BitmapData = new BitmapData(sw * DRAW_SCALE, sh * DRAW_SCALE, false, 0x0);
		
		private var drawMatrix:Matrix = new Matrix(DRAW_SCALE, 0, 0, DRAW_SCALE, 0, 0);
		private var drawColor:ColorTransform = new ColorTransform(0.1, 0.1, 0.1);
		
		private var sparkle:BitmapData = new BitmapData(sw/4, sh/4, true, 0x0);
		private var sparkleDrawMatrix:Matrix = new Matrix(0.25, 0, 0, 0.25, 0, 0);
		
		public static var mx:uint = 0;
		public static var my:uint = 0;
		
		public static var frameCount:uint = 0;
		
		private var display:Bitmap;
		private var pm:ParticleManager;
		private var prevMouse:Point = new Point();
		private var fluidImage:Bitmap;
		private var controls:Controls;
		
		public function Main():void 
		{
			if (stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		private function init(e:Event = null):void 
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			initStage();
			
			fSolver.rgb = true;
			fSolver.fadeSpeed = .007;
			fSolver.deltaT = .5;
			fSolver.viscosity = .00015;
			fSolver.vorticityConfinement = false;
			
			var dw:Number = sw / fluid.width;
			var dh:Number = sh / fluid.height;
			var s:Number = dw > dh ? dw : dh;
			
			pm = new ParticleManager();
			
			fluidImage = new Bitmap( fluid, 'never', true );
			var fadeImage:Bitmap = new Bitmap(fade, 'never', true);
			display = new Bitmap( screen, 'never', true );
			
			//fluidImage.scaleX = fluidImage.scaleY = s;
			fluidImage.width = sw;
			fluidImage.height = sh;
			fadeImage.scaleX = fadeImage.scaleY = 1 / DRAW_SCALE;
			
			display.blendMode = fadeImage.blendMode = BlendMode.ADD;
			
			fluidImage.y = display.y = fadeImage.y = 40;
			
			var b:Bitmap = new Bitmap(sparkle);
			b.y = 40;
			b.smoothing = true;
			b.blendMode = BlendMode.ADD;
			b.scaleX = b.scaleY = 4;
			
			controls = new Controls();
			
			addChild(fluidImage);
			addChild(fadeImage);
			addChild(display);
			addChild(b);
			addChild(controls);
			
			addEventListener(Event.ENTER_FRAME, render);
			stage.addEventListener(MouseEvent.MOUSE_MOVE, onMouseMove);
			stage.addEventListener(MouseEvent.CLICK, restart);
		}
		
		private function restart(e:MouseEvent):void
		{
			if (mouseY <= 40) return;
			
			drawMode = ++drawMode % 4;
			
			fade.fillRect(fade.rect, 0x0);
			screen.fillRect(screen.rect, 0x0);
			sparkle.fillRect(sparkle.rect, 0x0);
			fluid.fillRect(fluid.rect, 0x0);
			
			if(drawMode == 0) {
				drawParticles = true;
				drawFluid = true;
				drawLines = false;
			} else if (drawMode == 1) {
				drawParticles = false;
				drawFluid = true;
				drawLines = true;
			} else if (drawMode == 2) {
				drawParticles = false;
				drawFluid = true;
				drawLines = false;
			} else if(drawMode == 3){
				drawParticles = true;
				drawFluid = false;
				drawLines = false;
			}
		}
		
		private function onMouseMove(e:MouseEvent):void 
		{
			if (mouseY <= 40) return;
			
			handleForce(mouseX, mouseY - 40);
		}
		
		private function handleForce(x:Number, y:Number):void
		{
			const NormX:Number = x * isw;
			const NormY:Number = y * ish;
			const VelX:Number = (x - prevMouse.x) * isw;
			const VelY:Number = (y - prevMouse.y) * ish;

			addForce(NormX, NormY, VelX, VelY);
			
			prevMouse.x = x;
			prevMouse.y = y;
		}
		
		private function render(e:Event):void 
		{			
			//var t:int = getTimer();
			fSolver.update();
			//trace(getTimer() - t);
			
			if (drawFluid) {
				drawFluidBitmap();
			}
			
			if (drawParticles || drawLines) {
				//var t:int = getTimer();
				drawParticlesBitmap();
				//trace(getTimer() - t);
			}
			
			if(mx != 0 || my != 0) {
				moveScreen();
			}
			
			frameCount = ++frameCount % 0xFFFFFFFF;
		}
		
		public function drawFluidBitmap():void
		{
			const d:int = 0xFF * 1;
			const fw:int = fSolver.width;
			const tw:int = fw - 1;
			const th:int = fSolver.height - 1;
			
			var i:int, j:int, fi:int;
			var index:int = 0;
			
			for(j = 1; j < th; ++j) {
				for(i = 1; i < tw; ++i) {
					fi = int(i + fw * j);
					fluidBuffer[ index++ ] = ((fSolver.r[fi] * d) << 16) | ((fSolver.g[fi] * d) << 8) | (fSolver.b[fi] * d);
				}
			}
			fluid.lock();
			fluid.setVector( fluid.rect, fluidBuffer );
			fluid.applyFilter( fluid, fluid.rect, origin, blur );
			fluid.unlock(fluid.rect);
		}
		
		public function drawParticlesBitmap():void
		{
			screen.lock();
			if(drawLines) {
				screen.colorTransform( screen.rect, fade2alpha );
			} else {
				screen.fillRect(screen.rect, 0x0);
			}
			
			pm.update(screen, drawLines);
			screen.unlock( screen.rect );
			
			sparkle.lock();
			sparkle.fillRect(sparkle.rect, 0x0);
			sparkle.draw(screen, sparkleDrawMatrix);
			sparkle.unlock(sparkle.rect);
			
			if (frameCount & 1) {
				fade.lock();
				fade.draw(screen, drawMatrix, drawColor, BlendMode.ADD);
				fade.applyFilter(fade, fade.rect, origin, blur);
				fade.unlock(fade.rect);
			}
		}
		
		public function moveScreen():void
		{
			if (mx > 0 && frameCount & 1) fSolver.shiftLeft();
			if (my > 0 && frameCount & 1) fSolver.shiftUp();
			
			screen.scroll(-mx, -my);
			screen.fillRect(new Rectangle(sw - 4, 0, 4, sh), 0x000000);
		}
		
		public function addForce(x:Number, y:Number, dx:Number, dy:Number):void
		{
			const speed:Number = dx * dx  + dy * dy * aspectRatio2;    // balance the x and y components of speed with the screen aspect ratio
			
			if(speed > 0) {
				if (x < 0) x = 0;
				else if (x > 1) x = 1;
				if (y < 0) y = 0;
				else if (y > 1) y = 1;
				
				const colorMult:Number = 30;
				const velocityMult:Number = 20.0;
				
				const index:int = fSolver.getIndexForNormalizedPosition(x, y);
				
				const hue:Number = ((x + y) * 180 + frameCount) % 360;
				
				const rgb:Object = ColorUtils.HSB2GRB(hue);
				
				fSolver.rOld[index]  += rgb.r * colorMult;
				fSolver.gOld[index]  += rgb.g * colorMult;
				fSolver.bOld[index]  += rgb.b * colorMult;

				if(drawParticles || drawLines) pm.addParticles(x * sw, y * sh, 10);
				
				fSolver.uOld[index] += dx * velocityMult;
				fSolver.vOld[index] += dy * velocityMult;
			}
		}
		
		public static function setScroll(sx:Boolean = false, sy:Boolean = false):void
		{
			if (sx) {
				mx = 4;
			} else {
				mx = 0;
			}
			
			if (sy) {
				my = 4;
			} else {
				my = 0;
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