package ru.inspirit.utils 
{
	/***********************************************************************
 
	* This is a class for solving real-time fluid dynamics simulations based on Navier-Stokes equations 
	* and code from Jos Stam's paper "Real-Time Fluid Dynamics for Games" http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
	* Other useful resources and implementations I looked at while building this lib: 
	* Mike Ash (C) - http://mikeash.com/?page=pyblog/fluid-simulation-for-dummies.html
	* Alexander McKenzie (Java) - http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm
	* Pierluigi Pesenti (AS3 port of Alexander's) - http://blog.oaxoa.com/2008/01/21/actionscript-3-fluids-simulation/
	* Gustav Taxen (C) - http://www.nada.kth.se/~gustavt/fluids/
	* Dave Wallin (C++) - http://nuigroup.com/touchlib/ (uses portions from Gustav's)
	
	/***********************************************************************
 
	Copyright (c) 2008, 2009, Memo Akten, www.memo.tv
	*** The Mega Super Awesome Visuals Company ***

	/**
	 * AS3 Alchemy Port
	 * @author Eugene Zatepyakin
	 * @link http://blog.inspirit.ru/?p=339
	 * @link http://code.google.com/p/in-spirit/source/browse/#svn/trunk/projects/FluidSolverHD
	 */
	 	
	import cmodule.fluidsolver.CLibInit;
	import com.joa_ebert.apparat.memory.Memory;
	import flash.display.Bitmap;
	import flash.display.BitmapData;
	import flash.filters.BlurFilter;
	import flash.geom.Point;
	import flash.utils.ByteArray;
	
	public final class FluidSolverHD extends Bitmap 
	{
		protected const BLUR:BlurFilter = new BlurFilter(2, 2, 2);
		protected const ORIGIN:Point = new Point();
		
		protected const lib:Object = (new CLibInit()).init();
		
		protected var alchemyRAM:ByteArray;
		
		protected var imagePos:int;
		protected var particlesPos:int;
		protected var particlesPoolPos:int;
		protected var rOldPos:int;
		protected var gOldPos:int;
		protected var bOldPos:int;
		protected var uOldPos:int;
		protected var vOldPos:int;
		protected var particlesNumPos:int;
		protected var particlesDataPos:int;
		
		protected var width2:int;
		protected var height2:int;
		protected var screenW:int;
		protected var screenH:int;
		
		protected var _drawMode:int = 0;
		
		public var fluidImage:BitmapData;
		public var particlesImage:BitmapData;
		
		protected var clearParticles:ByteArray;
		protected var particlesArea:int;
		
		public function FluidSolverHD(width:int, height:int, screenW:int = 800, screenH:int = 600) 
		{
			fluidImage = new BitmapData(width, height, false, 0x00);
			particlesImage = new BitmapData(screenW, screenH, true, 0);
			
			super(fluidImage, 'never', true);
			
			width2 = width + 2;
			height2 = height + 2;
			this.screenW = screenW;
			this.screenH = screenH;
			
			clearParticles = particlesImage.getPixels(particlesImage.rect);
			clearParticles.position = 0;
			
			var ns : Namespace = new Namespace( "cmodule.fluidsolver" );
			alchemyRAM = (ns::gstate).ds;
			
			imagePos = lib.setupSolver(width, height, screenW, screenH);
			particlesPos = lib.getParticlesPointer();
			particlesPoolPos = lib.getParticlesPoolPointer();
			rOldPos = lib.getROldPointer();
			gOldPos = lib.getGOldPointer();
			bOldPos = lib.getBOldPointer();
			uOldPos = lib.getUOldPointer();
			vOldPos = lib.getVOldPointer();
			particlesNumPos = lib.getParticlesCountPointer();
			particlesDataPos = lib.getParticlesDataPointer();
		}

		public function update():void
		{
			/*if(_drawMode > 0){
				alchemyRAM.position = particlesPos;
				alchemyRAM.writeBytes(clearParticles);
			}*/
			
			lib.updateSolver();
			
			alchemyRAM.position = imagePos;
			fluidImage.lock();
			fluidImage.setPixels(fluidImage.rect, alchemyRAM);
			fluidImage.applyFilter(fluidImage, fluidImage.rect, ORIGIN, BLUR);
			fluidImage.unlock(fluidImage.rect);
			
			if(_drawMode > 0){
				drawParticles();
				//alchemyRAM.position = particlesPos;
				//particlesImage.lock();
				//particlesImage.setPixels(particlesImage.rect, alchemyRAM);
				//particlesImage.unlock(particlesImage.rect);
			}
		}
		
		public function drawParticles():void
		{
			var pos:int = Memory.readInt(particlesDataPos);
			var pn:int = Memory.readInt(particlesNumPos);
			var step:int = 6 << 3;
			var aa:int, xp:int, yp:int, vx:int, vy:int, cc:uint;
			var x:int, y:int, dx:int, dy:int, xinc:int, yinc:int, cumul:int, i:int;			
			
			particlesImage.lock();
			particlesImage.fillRect(particlesImage.rect, 0);
			
			while( --pn > -1 )
			{
				aa = int( Memory.readDouble(pos + 0)*0xFF + 0.5 );
				xp = int( Memory.readDouble(pos + 8) + 0.5 );
				yp = int( Memory.readDouble(pos + 16) + 0.5 );
				vx = int( Memory.readDouble(pos + 24) + 0.5 );
				vy = int( Memory.readDouble(pos + 32) + 0.5 );
				cc = (aa<<24) | (aa<<16) | (aa<<8) | aa;
				
				// inlined line drawing method
				
				x = xp - vx;
				y = yp - vy;		
				dx = xp - x;
				dy = yp - y;
				xinc = ( dx > 0 ) ? 1 : -1;
				yinc = ( dy > 0 ) ? 1 : -1;			
				
				dx = (dx ^ (dx >> 31)) - (dx >> 31);
				dy = (dy ^ (dy >> 31)) - (dy >> 31);
				
				particlesImage.setPixel32(x, y, cc);
				
				if ( dx > dy ) {
					cumul = dx >> 1 ;
					for ( i = 1;i <= dx; ++i ) 
					{
						x += xinc;
						cumul += dy;
						if (cumul >= dx) {
							cumul -= dx;
							y += yinc;
						}
						particlesImage.setPixel32(x, y, cc);
					}
				} else {
					cumul = dy >> 1;
					for ( i = 1; i <= dy; ++i ) 
					{
						y += yinc;
						cumul += dx;
						if ( cumul >= dy ) {
							cumul -= dy;
							x += xinc;
						}
						particlesImage.setPixel32(x, y, cc);
					}
				}
				
				pos += step;
			}
			
			particlesImage.unlock(particlesImage.rect);
		}

		public function addForce(tx:Number, ty:Number, dx:Number, dy:Number, rgb:Object, addColor:Boolean = true, addForce:Boolean = false, colorMult:Number = 50, velocityMult:Number = 30):void
		{
			var nx:int = int(tx * width2);
			var ny:int = int(ty * height2);
			if(nx < 1) nx=1; else if(nx > width) nx = width;
			if(ny < 1) ny=1; else if(ny > height) ny = height;
			
			var index:int = (nx + width2 * ny) << 3;
			
			if(addColor){
				Memory.writeDouble( Memory.readDouble(rOldPos + index) + rgb.r * colorMult, rOldPos + index);
				Memory.writeDouble( Memory.readDouble(gOldPos + index) + rgb.g * colorMult, gOldPos + index);
				Memory.writeDouble( Memory.readDouble(bOldPos + index) + rgb.b * colorMult, bOldPos + index);
			}
			if(addForce){
				Memory.writeDouble( Memory.readDouble(uOldPos + index) + dx * velocityMult, uOldPos + index);
				Memory.writeDouble( Memory.readDouble(vOldPos + index) + dy * velocityMult, vOldPos + index);
			}
			
			// Add some particles
			if(_drawMode > 0) addParticles(tx * particlesImage.width, ty * particlesImage.height, particlesPoolPos);
			
			/*alchemyRAM.position = rOldPos + index;
			var old:Number = alchemyRAM.readDouble();
			alchemyRAM.position = rOldPos + index;
			alchemyRAM.writeDouble(old + rgb.r * colorMult);
			
			alchemyRAM.position = gOldPos + index;
			old = alchemyRAM.readDouble();
			alchemyRAM.position = gOldPos + index;
			alchemyRAM.writeDouble(old + rgb.g * colorMult);
			
			alchemyRAM.position = bOldPos + index;
			old = alchemyRAM.readDouble();
			alchemyRAM.position = bOldPos + index;
			alchemyRAM.writeDouble(old + rgb.b * colorMult);
			
			alchemyRAM.position = uOldPos + index;
			old = alchemyRAM.readDouble();
			alchemyRAM.position = uOldPos + index;
			alchemyRAM.writeDouble(old + dx * velocityMult);
			
			alchemyRAM.position = vOldPos + index;
			old = alchemyRAM.readDouble();
			alchemyRAM.position = vOldPos + index;
			alchemyRAM.writeDouble(old + dy * velocityMult);*/
		}
		
		public function addParticles(tx:Number, ty:Number, pos:int) : void 
		{
			Memory.writeDouble(tx, pos);
			Memory.writeDouble(ty, pos + 8);
		}
		
		public function setWrap(wx:Boolean, wy:Boolean):void
		{
			lib.setWrap(wx ? 1 : 0, wy ? 1 : 0);
		}
		
		public function set colorDiffusion(value:Number):void
		{
			lib.setcolorDiffusion(value);
		}
		
		public function set solverIterations(value:int):void
		{
			lib.setsolverIterations(value);
		}
		
		public function set viscosity(value:Number):void
		{
			lib.setViscosity(value);
		}

		public function set fadeSpeed(value:Number):void
		{
			lib.setFadeSpeed(value);
		}
		
		public function set deltaT(value:Number):void
		{
			lib.setDeltaT(value);
		}
		
		public function set drawMode(value:int):void
		{
			lib.setDrawMode(value);
			_drawMode = value;
			visible = value < 3;
			if(_drawMode == 0){
				lib.clearParticles();
			}
		}
		
		public function set vorticityConfinement(value:Boolean):void
		{
			lib.setVorticityConfinement( value ? 1 : 0 );
		}
		
		public function get drawMode():int
		{
			return _drawMode;
		}
		
		public function get particlesNumber():int
		{
			return Memory.readInt(particlesNumPos);
		}
	}
}
