#pragma once
namespace Engine
{
	namespace Readers
	{
		void b64_decode_mio (char *dest, char *src);
		int b64_encode (char *dest,	const unsigned char *src,	int len);
		int b64_decode (char *dest,	const char *src);
	}
}
