#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "hdf.h"
#include "mfhdf.h"

static char *h4_string(const char *value, int32_t length)
{
    char *result = malloc((size_t)length + 1);
    if (result == NULL)
        return NULL;

    memcpy(result, value, (size_t)length);
    result[length] = '\0';
    return result;
}

static size_t h4_string_length(const char *value, size_t maximum)
{
    size_t length = 0;
    while (length < maximum && value[length] != '\0')
        ++length;
    return length;
}

int32_t ldas_h4_hopen(const char *path, int32_t path_length, int32_t access, int32_t ndds)
{
    char *c_path = h4_string(path, path_length);
    int32_t file_id = FAIL;

    if (c_path != NULL) {
        file_id = Hopen(c_path, (intn)access, (int16)ndds);
        free(c_path);
    }
    return file_id;
}

int32_t ldas_h4_vfstart(int32_t file_id) { return Vinitialize(file_id); }
int32_t ldas_h4_vsfatch(int32_t file_id, int32_t reference, const char *access, int32_t access_length)
{
    char *c_access = h4_string(access, access_length);
    int32_t vdata_id = FAIL;

    if (c_access != NULL) {
        vdata_id = VSattach(file_id, reference, c_access);
        free(c_access);
    }
    return vdata_id;
}
int32_t ldas_h4_vsqfnelt(int32_t vdata_id, int32_t *records) { *records = VSelts(vdata_id); return *records < 0 ? FAIL : SUCCEED; }
int32_t ldas_h4_vsfseek(int32_t vdata_id, int32_t position) { return VSseek(vdata_id, position); }
int32_t ldas_h4_vsfsfld(int32_t vdata_id, const char *fields, int32_t fields_length)
{
    char *c_fields = h4_string(fields, fields_length);
    int32_t status = FAIL;

    if (c_fields != NULL) {
        status = VSsetfields(vdata_id, c_fields);
        free(c_fields);
    }
    return status;
}
int32_t ldas_h4_vsfread(int32_t vdata_id, void *data, int32_t records, int32_t interlace)
{
    return VSread(vdata_id, data, records, interlace);
}
int32_t ldas_h4_vsfdtch(int32_t vdata_id) { return VSdetach(vdata_id); }
int32_t ldas_h4_vfend(int32_t file_id) { return Vfinish(file_id); }
int32_t ldas_h4_hclose(int32_t file_id) { return Hclose(file_id); }

int32_t ldas_h4_sfstart(const char *path, int32_t path_length, int32_t access)
{
    char *c_path = h4_string(path, path_length);
    int32_t sd_id = FAIL;

    if (c_path != NULL) {
        sd_id = SDstart(c_path, access);
        free(c_path);
    }
    return sd_id;
}
int32_t ldas_h4_sfn2index(int32_t sd_id, const char *name, int32_t name_length)
{
    char *c_name = h4_string(name, name_length);
    int32_t index = FAIL;

    if (c_name != NULL) {
        index = SDnametoindex(sd_id, c_name);
        free(c_name);
    }
    return index;
}
int32_t ldas_h4_sfselect(int32_t sd_id, int32_t index) { return SDselect(sd_id, index); }
int32_t ldas_h4_sfginfo(int32_t sds_id, char *name, int32_t name_length, int32_t *rank, int32_t *dimsizes, int32_t *data_type, int32_t *attributes)
{
    char *c_name = calloc((size_t)name_length + 1, 1);
    int32 hdf_dimsizes[H4_MAX_VAR_DIMS];
    int32_t status = FAIL;

    if (c_name != NULL) {
        status = SDgetinfo(sds_id, c_name, rank, hdf_dimsizes, data_type, attributes);
        if (status == SUCCEED) {
            memcpy(name, c_name, h4_string_length(c_name, (size_t)name_length));
            dimsizes[0] = hdf_dimsizes[0];
            dimsizes[1] = hdf_dimsizes[1];
        }
        free(c_name);
    }
    return status;
}
int32_t ldas_h4_sfrdata(int32_t sds_id, const int32_t *start, const int32_t *stride, const int32_t *edge, void *data)
{
    return SDreaddata(sds_id, (int32_t *)start, (int32_t *)stride, (int32_t *)edge, data);
}
int32_t ldas_h4_sfendacc(int32_t sds_id) { return SDendaccess(sds_id); }
int32_t ldas_h4_sfend(int32_t sd_id) { return SDend(sd_id); }
