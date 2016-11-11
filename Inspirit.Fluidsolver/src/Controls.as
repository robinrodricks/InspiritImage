package  
{
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
		
		private var p:Panel;
		
		public function Controls() 
		{
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
			p.width = 800;
			p.height = 40;
			
			var lb:Label = new Label(p, 10, 5);
			lb.name = 'fps_txt';
			
			lb = new Label(this, 40, 5, 'CLICK TO SWITCH BETWEEN RENDER MODES');
			lb.x = 10;
			lb.y = Main.sh + 15;
			
			var sl:HUISlider = new HUISlider(p, 90, 5, 'DELTA ', onDeltaChange);
			sl.setSliderParams(.2, 2, .5);
			sl.width = 180;
			
			sl = new HUISlider(p, 90, 18, 'VISCOS', onVicosityChange);
			sl.setSliderParams(.0001, .002, .00015);
			sl.labelPrecision = 4;
			sl.width = 180;
			
			sl = new HUISlider(p, 265, 5, 'FADE', onFadeChange);
			sl.setSliderParams(.001, .01, .007);
			sl.labelPrecision = 3;
			sl.width = 180;
			
			sl = new HUISlider(p, 265, 18, 'PREC', onIterChange);
			sl.setSliderParams(1, 10, 10);
			sl.labelPrecision = 0;
			sl.width = 180;
			
			lb = new Label(p, 440, 5, 'COLOR DIFFUSION');
			sl = new HUISlider(p, 433, 18, '', onColorDiffChange);
			sl.setSliderParams(0, .001, 0);
			sl.labelPrecision = 4;
			sl.width = 150;
			
			lb = new Label(p, 580, 5, 'VORTICITY CONF');
			var chk:CheckBox = new CheckBox(p, 580, 22, 'ENABLE', onVorticityChange);
			chk.name = 'vort';
			
			chk = new CheckBox(p, 675, 9, 'WRAP-X', onWrapChange);
			chk.name = 'wx';
			chk = new CheckBox(p, 675, 22, 'WRAP-Y', onWrapChange);
			chk.name = 'wy';
			
			chk = new CheckBox(p, 730, 9, 'SCROLL-X', onScrollChange);
			chk.name = 'sx';
			chk = new CheckBox(p, 730, 22, 'SCROLL-Y', onScrollChange);
			chk.name = 'sy';
			
			addEventListener(Event.ENTER_FRAME, countFrameTime);
		}
		
		private function onVorticityChange(e:Event):void
		{
			Main.fSolver.vorticityConfinement = CheckBox(p.getChildByName('vort')).selected;
		}
		
		private function onScrollChange(e:Event):void
		{
			var sx:Boolean = CheckBox(p.getChildByName('sx')).selected;
			var sy:Boolean = CheckBox(p.getChildByName('sy')).selected;
			Main.setScroll(sx, sy);
		}
		
		private function onWrapChange(e:Event):void
		{
			var wx:Boolean = CheckBox(p.getChildByName('wx')).selected;
			var wy:Boolean = CheckBox(p.getChildByName('wy')).selected;
			Main.fSolver.setWrap(wx, wy);
		}
		
		private function onColorDiffChange(e:Event):void
		{
			Main.fSolver.colorDiffusion = HUISlider(e.target).value;
		}
		
		private function onIterChange(e:Event):void
		{
			Main.fSolver.solverIterations = HUISlider(e.target).value;
		}
		
		private function onFadeChange(e:Event):void
		{
			Main.fSolver.fadeSpeed = HUISlider(e.target).value;
		}
		
		private function onVicosityChange(e:Event):void
		{
			Main.fSolver.viscosity = HUISlider(e.target).value;
		}
		
		private function onDeltaChange(e:Event):void
		{
			Main.fSolver.deltaT = HUISlider(e.target).value;
		}
		
		private function countFrameTime(e:Event = null):void
		{
			_timer = getTimer();
			if( _timer - 1000 >= _ms_prev )
			{
				_ms_prev = _timer;
				
				Label(p.getChildByName('fps_txt')).text = 'FPS: ' + _fps + ' / ' + stage.frameRate + '\nMS:';
				
				_fps = 0;
			}
			
			_fps ++;
			Label(p.getChildByName('fps_txt')).text = Label(p.getChildByName('fps_txt')).text.split('MS:')[0] + 'MS: ' + (_timer - _ms);
			_ms = _timer;
		}
		
	}
	
}