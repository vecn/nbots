#include <jni.h>

static const char* get_val(int status);

jobject jModelStatus_create(JNIEnv *env, int status)
{
	const char *str_class = "nb/geometricBot/ModelStatus";
	jclass class = (*env)->FindClass(env, str_class);
	const char *val = get_val(status);
	jfieldID field_id =
		(*env)->GetStaticFieldID(env, class, val,
					 "Lnb/geometricBot/ModelStatus;");
	jobject instance =
		(*env)->GetStaticObjectField(env, class, field_id);
	return instance;
}

static const char* get_val(int status)
{
	const char* val;
	switch (status) {
	case 0:
		val = "OK";
		break;
	case 1:
		val = "ZERO_VERTICES";
		break;
	case 2:
		val = "ZERO_EDGES";
		break;
	case 3:
		val = "REPEATED_VERTICES";
		break;
	case 4:
		val = "INCOHERENT_EDGES";
		break;
	case 5:
		val = "REPEATED_EDGES";
		break;
	case 6:
		val = "INTERSECTED_EDGES";
		break;
	case 7:
		val = "VTX_INTERSECTING_EDGE";
		break;
	case 8:
		val = "UNCLOSED_BOUNDARY";
		break;
	default:
		val = "OK";
	}
	return val;
}
