package ru.inspirit.surf 
{
	import apparat.memory.Memory;

	import cmodule.assurf_clib.CLibInit;

	import ru.inspirit.surf.ar.ARCamera3D;

	import flash.display.BitmapData;
	import flash.display.Graphics;
	import flash.display.Shape;
	import flash.filters.ConvolutionFilter;
	import flash.geom.Matrix;
	import flash.geom.Point;
	import flash.geom.Rectangle;
	import flash.system.Capabilities;
	import flash.utils.ByteArray;

	/**
	 * @author Eugene Zatepyakin
	 */
	public final class ASSURF 
	{		
		protected static const ASSURF_LIB:Object = (new CLibInit()).init();
		protected static const ns:Namespace = new Namespace( "cmodule.assurf_clib" );
		protected static const alchemyRAM:ByteArray = (ns::gstate).ds;
		
		public static const GET_MATCHES:int = 0;
		public static const GET_HOMOGRAPHY:int = 1;
		public static const GET_3DPOSE:int = 2;
		
		public static const DETECT_PRECISION_HIGH:int = 0;
		public static const DETECT_PRECISION_MEDIUM:int = 1;
		public static const DETECT_PRECISION_LOW:int = 2;
		
		protected static const GAUSSIAN_3x3:ConvolutionFilter = new ConvolutionFilter(
																				3,3,[ 1,2,1,
																				2,4,2,
																				1,2,1], 16
																				);
		protected static const ORIGIN:Point = new Point();
		protected static const DEF_IMAGE_PROCESSOR:ImageProcessor = new ImageProcessor();
		
		protected var maxScreenPoints:int;
		
		protected var imgProcessor:ImageProcessor = DEF_IMAGE_PROCESSOR;
		protected var buffer:BitmapData;
		protected var rect:Rectangle;
		
		protected var rect2:Rectangle;
		protected var rect3:Rectangle;
		
		protected var buffer2:BitmapData;
		protected var scale2_mtx:Matrix = new Matrix(0.5, 0, 0, 0.5);
		protected var scale4_mtx:Matrix = new Matrix(0.25, 0, 0, 0.25);
		protected var img2Ptr:int;
		protected var blur2Ptr:int;
		protected var img3Ptr:int;
		protected var blur3Ptr:int;
		
		protected var useMask:int = 0;
		protected var mask_bmp:BitmapData;
		protected var mask_mtx:Matrix;
		protected const mask_sh:Shape = new Shape();
		protected const mask_gfx:Graphics = mask_sh.graphics;

		protected var testPtr:int;
		protected var imgPtr:int;
		protected var blurPtr:int;
		protected var maskPtr:int;
		protected var arCameraPtr:int;
		protected var refIndexesPtr:int;
		protected var supressNeighbPtr:int;
		
		protected var referenceMap:Vector.<ReferenceInfo>;
		
		public var autoDetectROI:Boolean = false;
		
		public function init(detectPrecision:int = DETECT_PRECISION_MEDIUM, maxScreenPoints:int = 500, maxRefPoints:int = 5000, maxRefObjects:int = 1):void
		{
			ASSURF_LIB.setupGlobalBuffers(maxRefPoints, maxRefObjects, maxScreenPoints, detectPrecision);
			
			this.maxScreenPoints = maxScreenPoints;
			
			if(referenceMap) clearRefObjects();
			referenceMap = new Vector.<ReferenceInfo>(maxRefObjects, true);
			updateDataPointers();
		}
		
		public function setup(width:int, height:int):void
		{
			buffer = new BitmapData(width, height, false, 0x00);
			mask_bmp = buffer.clone();
			rect = buffer.rect;
			imgProcessor.imageRect = rect;
			DEF_IMAGE_PROCESSOR.imageRect = rect;
			
			buffer2 = new BitmapData( width, height, false, 0x00 );
			rect2 = new Rectangle( 0, 0, width >> 1, height >> 1);
			rect3 = new Rectangle( 0, 0, width >> 2, height >> 2);
			
			mask_mtx = new Matrix(1.2, 0, 0, 1.2, -(width*0.5)*0.2, -(height*0.5)*0.2);
			
			var ptrs : Array = ASSURF_LIB.setupImagePyramid(width, height);
			
			imgPtr = ptrs[0];
			img2Ptr = ptrs[1];
			img3Ptr = ptrs[2];
			blurPtr = ptrs[3];
			blur2Ptr = ptrs[4];
			blur3Ptr = ptrs[5];
			maskPtr = ptrs[6];
		}
		
		public function detectSingleObject(bmp:BitmapData, objID:int = 0, options:int = GET_MATCHES):ReferenceInfo
		{
			imgProcessor.preProcess(bmp, buffer);
			
			var i:int;
			
			var data1:Vector.<uint> = buffer.getVector(rect);
			

			buffer2.draw( buffer, scale2_mtx, null, null, rect, true );
			var data2:Vector.<uint> = buffer2.getVector(rect2);
			
			buffer2.applyFilter(buffer2, rect2, ORIGIN, GAUSSIAN_3x3);
			var data2_b:Vector.<uint> = buffer2.getVector(rect2);
			
			
			buffer2.draw( buffer, scale4_mtx, null, null, rect, true );
			var data3:Vector.<uint> = buffer2.getVector(rect3);
			
			buffer2.applyFilter(buffer2, rect3, ORIGIN, GAUSSIAN_3x3);
			var data3_b:Vector.<uint> = buffer2.getVector(rect3);
			
			
			buffer.applyFilter(buffer, rect, ORIGIN, GAUSSIAN_3x3);
			var data1_b:Vector.<uint> = buffer.getVector(rect);
			
			
			var len1:int = data1.length;
			var len2:int = len1 >> 2;
			var len3:int = len2 >> 2;
			var addr_o1:int = imgPtr;
			var addr_b1:int = blurPtr;
			var addr_o2:int = img2Ptr;
			var addr_b2:int = blur2Ptr;
			var addr_o3:int = img3Ptr;
			var addr_b3:int = blur3Ptr;
			
			for(i = 0; i < len3; ++i, ++addr_o1, ++addr_o2, ++addr_o3, ++addr_b1, ++addr_b2, ++addr_b3)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
				//
				Memory.writeByte(data2[i] & 0xFF, addr_o2);
				Memory.writeByte(data2_b[i] & 0xFF, addr_b2);
				//
				Memory.writeByte(data3[i] & 0xFF, addr_o3);
				Memory.writeByte(data3_b[i] & 0xFF, addr_b3);
			}
			for(; i < len2; ++i, ++addr_o1, ++addr_o2, ++addr_b1, ++addr_b2)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
				//
				Memory.writeByte(data2[i] & 0xFF, addr_o2);
				Memory.writeByte(data2_b[i] & 0xFF, addr_b2);
			}
			for(; i < len1; ++i, ++addr_o1, ++addr_b1)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
			}
			
			if(useMask == 1)
			{
				data1 = mask_bmp.getVector( rect );
				len1 = data1.length;
				addr_o1 = maskPtr;
				for(i = 0; i < len1; ++i, ++addr_o1) 
				{
					Memory.writeByte(data1[i] & 0xFF, addr_o1);
				}
			}
			
			/*var num_m:int = */ASSURF_LIB.runTask(maxScreenPoints, useMask, objID, options);
			
			var ref:ReferenceInfo = referenceMap[objID];
			
			ref.updateInfo();
			if(ref.matchedPointsCount > 4 && (autoDetectROI || options > 0))
			{
				ref.updateHomography();
			}
			
			if(options == 2 && ref.matchedPointsCount > 4) 
			{
				ref.updateTransform();
			}
			
			if(autoDetectROI)
			{
				updateMask(Vector.<ReferenceInfo>([ref]));
			}
			
			return ref;
		}
		
		public function detectMultiObject(bmp:BitmapData, options:int = GET_MATCHES):Vector.<ReferenceInfo>
		{			
			imgProcessor.preProcess(bmp, buffer);
			
			var i:int;
			var data1:Vector.<uint> = buffer.getVector(rect);
			

			buffer2.draw( buffer, scale2_mtx, null, null, rect, true );
			var data2:Vector.<uint> = buffer2.getVector(rect2);
			
			buffer2.applyFilter(buffer2, rect2, ORIGIN, GAUSSIAN_3x3);
			var data2_b:Vector.<uint> = buffer2.getVector(rect2);
			
			
			buffer2.draw( buffer, scale4_mtx, null, null, rect, true );
			var data3:Vector.<uint> = buffer2.getVector(rect3);
			
			buffer2.applyFilter(buffer2, rect3, ORIGIN, GAUSSIAN_3x3);
			var data3_b:Vector.<uint> = buffer2.getVector(rect3);
			
			
			buffer.applyFilter(buffer, rect, ORIGIN, GAUSSIAN_3x3);
			var data1_b:Vector.<uint> = buffer.getVector(rect);
			
			
			var len1:int = data1.length;
			var len2:int = len1 >> 2;
			var len3:int = len2 >> 2;
			var addr_o1:int = imgPtr;
			var addr_b1:int = blurPtr;
			var addr_o2:int = img2Ptr;
			var addr_b2:int = blur2Ptr;
			var addr_o3:int = img3Ptr;
			var addr_b3:int = blur3Ptr;
			
			for(i = 0; i < len3; ++i, ++addr_o1, ++addr_o2, ++addr_o3, ++addr_b1, ++addr_b2, ++addr_b3)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
				//
				Memory.writeByte(data2[i] & 0xFF, addr_o2);
				Memory.writeByte(data2_b[i] & 0xFF, addr_b2);
				//
				Memory.writeByte(data3[i] & 0xFF, addr_o3);
				Memory.writeByte(data3_b[i] & 0xFF, addr_b3);
			}
			for(; i < len2; ++i, ++addr_o1, ++addr_o2, ++addr_b1, ++addr_b2)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
				//
				Memory.writeByte(data2[i] & 0xFF, addr_o2);
				Memory.writeByte(data2_b[i] & 0xFF, addr_b2);
			}
			for(; i < len1; ++i, ++addr_o1, ++addr_b1)
			{
				Memory.writeByte(data1[i] & 0xFF, addr_o1);
				Memory.writeByte(data1_b[i] & 0xFF, addr_b1);
			}
			
			if(useMask == 1)
			{
				data1 = mask_bmp.getVector( rect );
				len1 = data1.length;
				addr_o1 = maskPtr;
				for(i = 0; i < len1; ++i, ++addr_o1) 
				{
					Memory.writeByte(data1[i] & 0xFF, addr_o1);
				}
			}
			
			var num_m:int = ASSURF_LIB.runTask2(maxScreenPoints, useMask, options);
			
			var res:Vector.<ReferenceInfo> = new Vector.<ReferenceInfo>(num_m, true);
			
			for(i = 0; i < num_m; ++i)
			{
				res[i] = getReferenceInfo(Memory.readInt(refIndexesPtr + (i << 2)), (autoDetectROI || options > 0), options == 2);
			}
			
			if(autoDetectROI)
			{
				updateMask(res);
			}
			
			return res;
		}
		
		public function getReferenceInfo(refID:int = 0, updateHomography:Boolean = false, updateTransform:Boolean = false):ReferenceInfo
		{
			var ref:ReferenceInfo = referenceMap[refID];
			
			ref.updateInfo();
			if(updateHomography && ref.matchedPointsCount > 4) 
			{
				ref.updateHomography();
			}
			
			if(updateTransform && ref.matchedPointsCount > 4) 
			{
				ref.updateTransform();
			}
			
			return ref;
		}
		
		public function addRefObject(bmp:BitmapData, scaleLevels:uint = 4, maxPointsPerLevel:uint = 1000, supressNeighbors:Boolean = false):int
		{
			var iw:int = bmp.width;
			var ih:int = bmp.height;
			
			var objPtrs:Array = ASSURF_LIB.createRefObject(iw, ih);
			var objID:int = objPtrs[0];
			
			var ds:Number = 1.0;
			var ds_inc:Number = Math.SQRT2;

			Memory.writeInt(supressNeighbors ? 1 : 0, supressNeighbPtr);			
			
			for(var j:int = 0; j < scaleLevels; ++j)
			{
				var img_lev1:BitmapData = new BitmapData(iw/ds, ih/ds, false, 0x00);
				img_lev1.draw(bmp, new Matrix(1/ds, 0, 0, 1/ds, 0, 0), null, null, null, true);
				
				DEF_IMAGE_PROCESSOR.imageRect = img_lev1.rect;
				DEF_IMAGE_PROCESSOR.preProcess(img_lev1, img_lev1);
				
				var data:Vector.<uint> = img_lev1.getVector(img_lev1.rect);
				
				img_lev1.applyFilter(img_lev1, img_lev1.rect, ORIGIN, GAUSSIAN_3x3);
				var data_b:Vector.<uint> = img_lev1.getVector(img_lev1.rect);
				
				var ptrs : Array = ASSURF_LIB.setupImageHolders(img_lev1.width, img_lev1.height);
				imgPtr = ptrs[0];
				blurPtr = ptrs[1];
				
				var len:int = data.length;
				var addr:int = imgPtr;
				var addr_b:int = blurPtr;
				
				for(var i:int = 0; i < len; ++i, ++addr, ++addr_b) 
				{
					Memory.writeByte(data[i] & 0xFF, addr);
					Memory.writeByte(data_b[i] & 0xFF, addr_b);
				}
				
				ASSURF_LIB.pushImageToRefObject(objID, img_lev1.width, img_lev1.height, maxPointsPerLevel, bmp.width/img_lev1.width);
				
				ds *= ds_inc;
			}
			
			var objInf:ReferenceInfo = new ReferenceInfo(iw, ih);
			
			objInf.id = objID;
			objInf.pointsCountPtr = objPtrs[1];
			objInf.matchedPointsCountPtr = objPtrs[2];
			objInf.transformErrorPtr = objPtrs[3];
			objInf.homographyPtr = objPtrs[4];
			objInf.transformPtr = objPtrs[5];
			
			objInf.updateInfo();
			
			referenceMap[objInf.id] = objInf;
			
			return objID;
		}
		
		public function clearRefObjects():void
		{
			ASSURF_LIB.clearReferenceObjects();
			
			var n:int = referenceMap.length;
			for(var i:int = 0; i < n; ++i)
			{
				referenceMap[i] = null;
			}
		}
		
		public function buildRefIndex():void
		{
			ASSURF_LIB.buildRefIndex();
		}
		
		public function setupARCamera(arCam:ARCamera3D):void
		{
			Memory.writeDouble(arCam.m02, arCameraPtr);
			Memory.writeDouble(arCam.m12, arCameraPtr + 8);
			Memory.writeDouble(arCam.m00, arCameraPtr + 16);
			Memory.writeDouble(arCam.m11, arCameraPtr + 24);
		}
		
		public function exportReferenceData(ba:ByteArray):void 
		{
			ASSURF_LIB.exportReferencesData(ba);
		}
		
		public function importReferenceData(ba:ByteArray):void
		{
			clearRefObjects();
			
			/*var detectPrec:int = */ba.readInt();
			var descrSize:int = ba.readInt();
			var numRefs : int = ba.readInt();
			var refPointsCount:int = 0;

			//trace('data:', detectPrec, descrSize, numRefs);

			for(var i : int = 0; i < numRefs; ++i)
			{
				ba.position = 12 + (16 + (40 + (descrSize << 3)) * refPointsCount ) * i;
				/*var ind:int = */ba.readInt();
				var iw:int = ba.readInt();
				var ih:int = ba.readInt();
				var ptCnt:int = ba.readInt();

				//trace('ref:', ind, iw, ih, ptCnt);
				
				var objPtrs:Array = ASSURF_LIB.createRefObject(iw, ih);
				var objID:int = objPtrs[0];
				
				ASSURF_LIB.pushDataToRefObject(objID, ptCnt, ba);
				
				var objInf:ReferenceInfo = new ReferenceInfo(iw, ih);
			
				objInf.id = objID;
				objInf.pointsCountPtr = objPtrs[1];
				objInf.matchedPointsCountPtr = objPtrs[2];
				objInf.transformErrorPtr = objPtrs[3];
				objInf.homographyPtr = objPtrs[4];
				objInf.transformPtr = objPtrs[5];
				
				objInf.updateInfo();
				
				referenceMap[objInf.id] = objInf;
				
				refPointsCount += ptCnt;
			}
		}
		
		public function set imageProcessor(ip:ImageProcessor):void
		{
			this.imgProcessor = ip || DEF_IMAGE_PROCESSOR;
		}
		
		public function get ROImask():BitmapData
		{
			return mask_bmp;
		}
		
		public function get useROI():Boolean 
		{
			return useMask == 1;
		}
		
		public function set useROI(val:Boolean):void 
		{
			useMask = val ? 1 : 0;
		}

		protected function updateDataPointers():void
		{
			var pps:Array = ASSURF_LIB.getDataPointers();
			testPtr = pps[0];
			blurPtr = pps[1];
			maskPtr = pps[2];
			imgPtr = pps[3];
			arCameraPtr = pps[4];
			refIndexesPtr = pps[5];
			supressNeighbPtr = pps[6];

			img2Ptr = pps[7];
			blur2Ptr = pps[8];
			img3Ptr = pps[9];
			blur3Ptr = pps[10];
		}
		
		public function debug():String
		{
			var str:String = 'FP: ' + Capabilities.version + (Capabilities.isDebugger ? ' Debug' : '');
			str += '\nreference points: ' + Memory.readInt(testPtr + ((maxScreenPoints+3) << 2));
			str += '\nscreen points: ' + Memory.readInt(testPtr + ((maxScreenPoints)<<2));
			str += '\ncurr/prev frame matches: ' + Memory.readInt(testPtr + ((maxScreenPoints+1)<<2)) + '/' + Memory.readInt(testPtr + ((maxScreenPoints+2)<<2));
			str += '\nfiltered matches: ' + Memory.readInt(testPtr + ((maxScreenPoints+4)<<2)) + ' skip describe: ' + Boolean(Memory.readInt(testPtr + ((maxScreenPoints+5)<<2)));
			/*
			var pn:int = Memory.readInt(testPtr + ((maxScreenPoints)<<2));
			var ncc_arr:Array = [];
			for(var i:int = 0; i < pn; ++i)
			{
				var res:int = Memory.readInt(testPtr + (i<<2)); 
				ncc_arr.push(res / 1000);
			}
			ncc_arr.sort( Array.NUMERIC );
			trace('from-to', ncc_arr[0], ncc_arr[pn-1], '|', ncc_arr);
			*/
			
			str += '\ntracker: ' + Memory.readInt(testPtr + (0<<2)) + ' NCC calls, ';
			str += Memory.readInt(testPtr + (1<<2)) + ' ncc matches; ';
			str += 'tracked ' + Memory.readInt(testPtr + (2<<2)) + ' features with LK. Saved ' + Memory.readInt(testPtr + (3<<2)) + ' tracks.';
			
			return str;
		}
		public function get debugMask():Shape
		{
			return mask_sh;
		}
		
		public function get referencePointsCount():int
		{
			return Memory.readInt(testPtr + ((maxScreenPoints+3) << 2));
		}
		
		/**
		 * Clear all allocated memory inside Alchemy
		 */
		public function freeMemory():void 
		{
			ASSURF_LIB.clearBuffers();
		}

		protected function updateMask(objs:Vector.<ReferenceInfo>):void
		{
			var n:int = objs.length;
			var ref:ReferenceInfo;
			var h:HomographyMatrix;
			var good:int = 0;
			
			mask_bmp.fillRect(rect, 0x00);
			useMask = 0;
			mask_gfx.clear();
			
			for(var i:int = 0; i < n; ++i)
			{
				ref = objs[i];
				if(ref.matchedPointsCount > 4)
				{
					mask_gfx.beginFill(0x0000FF);
					
					h = ref.homography;
					var m11:Number = h.m11;
					var m12:Number = h.m12;
					var m13:Number = h.m13;
					var m21:Number = h.m21;
					var m22:Number = h.m22;
					var m23:Number = h.m23;
					var m31:Number = h.m31;
					var m32:Number = h.m32;
					var m33:Number = h.m33;
					var px:Number;
					var py:Number;
					var x:Number = ref.pt0.x;
					var y:Number = ref.pt0.y;
					var Z:Number = 1.0 / (m31 * x + m32 * y + m33);
					var p0x:Number = px = (m11*x + m12*y + m13) * Z;
					var p0y:Number = py = (m21*x + m22*y + m23) * Z;
					
					mask_gfx.moveTo(px, py);
					
					x = ref.pt1.x;
					y = ref.pt1.y;
					Z = 1.0 / (m31 * x + m32 * y + m33);
					px = (m11*x + m12*y + m13) * Z;
					py = (m21*x + m22*y + m23) * Z;
					
					mask_gfx.lineTo(px, py);
					
					x = ref.pt2.x;
					y = ref.pt2.y;
					Z = 1.0 / (m31 * x + m32 * y + m33);
					px = (m11*x + m12*y + m13) * Z;
					py = (m21*x + m22*y + m23) * Z;
					
					mask_gfx.lineTo(px, py);
					
					x = ref.pt3.x;
					y = ref.pt3.y;
					Z = 1.0 / (m31 * x + m32 * y + m33);
					px = (m11*x + m12*y + m13) * Z;
					py = (m21*x + m22*y + m23) * Z;
					
					mask_gfx.lineTo(px, py);
					mask_gfx.lineTo(p0x, p0y);
					
					mask_gfx.endFill();
					
					good++;
				}
			}
			
			if(good > n>>1) 
			{
				useMask = 1;
				mask_bmp.draw(mask_sh, mask_mtx);
			}
		}
	}
}
