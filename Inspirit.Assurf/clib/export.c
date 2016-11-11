
/* read, write, seek, close a File inside Alchemy, File = ByteArray */
static int readba(void *cookie, char *dst, int size)
{	
	return AS3_ByteArray_readBytes(dst, (AS3_Val)cookie, size);
}

static int writeba(void *cookie, const char *src, int size)
{
	return AS3_ByteArray_writeBytes((AS3_Val)cookie, (char *)src, size);
}

static fpos_t seekba(void *cookie, fpos_t offs, int whence)
{
	return AS3_ByteArray_seek((AS3_Val)cookie, offs, whence);
}

static int closeba(void *cookie)
{
	AS3_Val zero = AS3_Int(0);
	AS3_SetS((AS3_Val)cookie, "position", zero);
	AS3_Release(zero);
	return 0;
}

AS3_Val exportReferencesData(void *data, AS3_Val args)
{
	AS3_Val src;
	FILE *input;
	const int n = 64*1024;
	int out_size = 0;
	char out[n];
	int *iptr;
	double *dptr;
	
	AS3_ArrayValue( args, "AS3ValType", &src );
	
	input = funopen((void *)src, readba, writeba, seekba, closeba);
	
	iptr = (int*)out;
	
	*iptr++ = DETECT_PRECISION;
	*iptr++ = DESCRIPTOR_SIZE;
	*iptr++ = referenceCount;
	out_size += 12;
	
	int i, j, k;
	for(i = 0; i < referenceCount; i++)
	{
		const RefObject *obj = &refObjectsMap[i];
		*iptr++ = obj->index;
		*iptr++ = obj->width;
		*iptr++ = obj->height;
		*iptr++ = obj->pointsCount;
		
		out_size += 16;
		
		for(j = 0; j < obj->pointsCount; j++)
		{
			const IPoint *ipt = obj->points + j;
			
			*iptr++ = ipt->index;
			*iptr++ = ipt->refIndex;
			*iptr++ = ipt->pos;
			*iptr++ = ipt->x;
			*iptr++ = ipt->y;
			*iptr++ = 2;/* = ipt->scale*/;
			
			out_size += 24;
			
			dptr = (double*)(out+out_size);
			
			*dptr++ = ipt->score;
			*dptr++ = ipt->orientation;
			
			out_size += 16;
			
			k = DESCRIPTOR_SIZE;
			const double *descr = ipt->descriptor;
			while( --k > -1 )
			{
				*dptr++ = *descr++;
			}
			
			out_size += DESCRIPTOR_SIZE * 8;
			
			if(out_size >= n - 572 - 64)
			{
				fwrite(out, 1, out_size, input);
				out_size = 0;
			}
			
			iptr = (int*)(out + out_size);
		}
	}
	
	fclose( input );
	
	return 0;
}