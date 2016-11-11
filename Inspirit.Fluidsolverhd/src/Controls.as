package  
{
	import ru.inspirit.utils.FluidSolverHD;
	import com.bit101.components.CheckBox;
	import com.bit101.components.HUISlider;
	import com.bit101.components.Label;
	import com.bit101.components.Panel;
	import com.bit101.components.Style;

	import flash.display.Sprite;
	import flash.events.Event;
	import flash.utils.getTimer;

	/**
	 * ...
	 * @author Eugene Zatepyakin
	 */
	public class Controls extends Sprite
	{
		
		private var _timer : uint;
		private var _fps : uint;
		private var _ms : uint;
		private var _ms_prev : uint;
		
		public var alg_t:uint = 0;
		public var alg_n:uint = 0;
		
		private var p:Panel;
		
		protected var fSolver:FluidSolverHD;
		
		public function Controls(fs:FluidSolverHD) 
		{
			fSolver = fs;
			if (stage) init();
			else addEventListener(Event.ADDED_TO_STAGE, init);
		}
		
		private function init(e:Event = null):void 
		{
			removeEventListener(Event.ADDED_TO_STAGE, init);
			
			Style.PANEL = 0x333333;
			Style.BUTTON_FACE = 0x333333;
			Style.LABEL_TEXT = 0xF6F6F6;
			
			p = new Panel(this);
			p.width = 960;
			p.height = 40;
			
			var lb:Label = new Label(p, 10, 5);
			lb.name = 'fps_txt';
			
			lb = new Label(p, 115, 5);
			lb.name = 'inf_txt';
			
			lb = new Label(this, 40, 5, 'CLICK TO SWITCH BETWEEN RENDER MODES');
			lb.x = 10;
			lb.y = 580 - 25;
			
			var sl:HUISlider = new HUISlider(p, 210, 5, 'DELTA ', onDeltaChange);
			sl.setSliderParams(.2, 2, .5);
			sl.width = 180;
			
			sl = new HUISlider(p, 210, 18, 'VISCOS', onVicosityChange);
			sl.setSliderParams(.0001, .002, .00015);
			sl.labelPrecision = 4;
			sl.width = 180;
			
			sl = new HUISlider(p, 385, 5, 'FADE', onFadeChange);
			sl.setSliderParams(.001, .01, .007);
			sl.labelPrecision = 3;
			sl.width = 180;
			
			sl = new HUISlider(p, 385, 18, 'ITER', onIterChange);
			sl.setSliderParams(1, 10, 10);
			sl.labelPrecision = 0;
			sl.width = 180;
			
			lb = new Label(p, 560, 5, 'COLOR DIFFUSION');
			sl = new HUISlider(p, 553, 18, '', onColorDiffChange);
			sl.setSliderParams(0, .001, 0);
			sl.labelPrecision = 4;
			sl.width = 150;
			
			lb = new Label(p, 700, 5, 'VORTICITY CONF');
			var chk:CheckBox = new CheckBox(p, 700, 22, 'ENABLE', onVorticityChange);
			chk.name = 'vort';
			
			chk = new CheckBox(p, 795, 9, 'WRAP-X', onWrapChange);
			chk.name = 'wx';
			chk = new CheckBox(p, 795, 22, 'WRAP-Y', onWrapChange);
			chk.name = 'wy';
			
			addEventListener(Event.ENTER_FRAME, countFrameTime);
		}
		
		private function onVorticityChange(e:Event):void
		{
			fSolver.vorticityConfinement = CheckBox(p.getChildByName('vort')).selected;
		}
		
		private function onWrapChange(e:Event):void
		{
			var wx:Boolean = CheckBox(p.getChildByName('wx')).selected;
			var wy:Boolean = CheckBox(p.getChildByName('wy')).selected;
			fSolver.setWrap(wx, wy);
		}
		
		private function onColorDiffChange(e:Event):void
		{
			fSolver.colorDiffusion = HUISlider(e.target).value;
		}
		
		private function onIterChange(e:Event):void
		{
			fSolver.solverIterations = HUISlider(e.target).value;
		}
		
		private function onFadeChange(e:Event):void
		{
			fSolver.fadeSpeed = HUISlider(e.target).value;
		}
		
		private function onVicosityChange(e:Event):void
		{
			fSolver.viscosity = HUISlider(e.target).value;
		}
		
		private function onDeltaChange(e:Event):void
		{
			fSolver.deltaT = HUISlider(e.target).value;
		}
		
		private function countFrameTime(e:Event = null):void
		{
			_timer = getTimer();
			if( _timer - 1000 >= _ms_prev )
			{
				_ms_prev = _timer;
				
				Label(p.getChildByName('fps_txt')).text = 'FPS: ' + _fps + ' / ' + stage.frameRate + '\nRENDER TIME: ' + int(alg_t / alg_n + 0.5) + 'ms';
				
				_fps = alg_n = alg_t = 0;
			}
			
			_fps ++;
			_ms = _timer;
			
			Label(p.getChildByName('inf_txt')).text = 'PARTICLES: ' + fSolver.particlesNumber;
		}
		
	}
	
}