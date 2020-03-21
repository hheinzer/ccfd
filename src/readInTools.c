/*
 * readInTools.c
 *
 * Created: Sat 21 Mar 2020 10:46:32 AM CET
 * Author : hhh
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "main.h"

typedef struct cmd_t cmd_t;
struct cmd_t {
	char *key;
	char *value;
	cmd_t *next, *prev;
};

cmd_t *firstCmd;

/*
 * Read .ini file and parse each line into a cmd_t object. All cmd_t
 * objects are connected to a list of commands starting with "firstCmd".
 */
void fillCmds(char iniFileName[STRLEN])
{
	FILE *iniFile;
	printf("Reading parameter file: '%s'\n", iniFileName);
	iniFile = fopen(iniFileName, "r");
	if (iniFile == NULL) {
		printf("Abort! Could not find specified file.\n");
		exit(1);
	}

	int i, j;
	char line[2 * STRLEN], tmp[STRLEN];
	cmd_t *currCmd = NULL;
	while (fgets(line, sizeof(line), iniFile)) {
		cmd_t *aCmd = malloc(sizeof(cmd_t));

		/* read key and make it lowercase */
		j = 0;
		for (i = 0; i < strlen(line); ++i) {
			if (isalnum(line[i])) {
				tmp[j++] = tolower(line[i]);
			} else if ((line[i] == '!') || (line[i] == '#') ||
				   (line[i] == '=')) {
				break;
			}
		}
		tmp[j] = '\0';

		/* continue and read value if key was read, else free aCmd */
		if (strlen(tmp) > 0) {
			/* store the key */
			aCmd->key = malloc(strlen(tmp) + 1);
			strcpy(aCmd->key, tmp);
			aCmd->next = NULL;

			/* read value */
			j = 0;
			for (; i < strlen(line); ++i) {
				if (isalnum(line[i]) || (line[i] == '_') ||
				    (line[i] == '.') || (line[i] == ',') ||
				    (line[i] == '-')) {
					tmp[j++] = line[i];
				} else if ((line[i] == '!') || (line[i] == '#')) {
					break;
				}
			}
			tmp[j] = '\0';

			/* store value */
			aCmd->value = malloc(strlen(tmp) + 1);
			strcpy(aCmd->value, tmp);

			/* complete the chain links */
			if (firstCmd == NULL) {
				currCmd = firstCmd = aCmd;
				currCmd->prev = NULL;
			} else {
				currCmd->next = aCmd;
				currCmd->next->prev = currCmd;
				currCmd = currCmd->next;
			}

		} else {
			free(aCmd);
		}
	}
}

/*
 * Delete a single node of the command list.
 */
void deleteCmd(cmd_t *aCmd)
{
	/* remove the key and the value from memory */
	free(aCmd->key);
	free(aCmd->value);

	/* link up the previous command to the next, if it exists. Otherwise
	 * we are at the beginning of the list and need to move the first
	 * command one over. */
	if (aCmd->prev) {
		aCmd->prev->next = aCmd->next;
	} else {
		firstCmd = aCmd->next;
	}

	/* link up the next command to the previous command, if it exists */
	if (aCmd->next) {
		aCmd->next->prev = aCmd->prev;
	}
}

/*
 * Find the provided key in the list of commands, and return the address of
 * the corresponding value string. Return NULL if key was not found.
 */
char *findCmd(const char *key, char defMsg[8], const char *proposal)
{
	/* create a lowercase copy of the key */
	char *keyLower = malloc(strlen(key) + 1);
	strcpy(keyLower, key);
	for (int i = 0; i < strlen(keyLower); ++i) {
		keyLower[i] = tolower(keyLower[i]);
	}

	/* search for the first matching key in the list */
	char *value = NULL;
	cmd_t *currCmd;
	for (currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		if (!strcmp(keyLower, currCmd->key)) {
			/* the keys are the same, create enough space for the
			 * value to be stored and break out of the loop */
			value = malloc(strlen(currCmd->value) + 1);
			strcpy(value, currCmd->value);
			strcpy(defMsg, "*CUSTOM");
			deleteCmd(currCmd);
			break;
		}
	}

	/* if a proposal is not given, then the key is not optional and is
	 * missing in the .ini file */
	if (!value) {
		if (!proposal) {
			printf("Inifile missing necessary keyword: %s\n",
					keyLower);
			exit(1);
		} else {
			value = malloc(strlen(proposal) + 1);
			strcpy(value, proposal);
			strcpy(defMsg, "DEFAULT");
		}
	}
	free(keyLower);

	return value;
}

/*
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 */
char *getStr(const char *key, const char *proposal)
{
	char defMsg[8];
	char *value = findCmd(key, defMsg, proposal);
	printf("| %19s = %27s (%s)\n", key, value, defMsg);
	return value;
}

/*
 * Count all occourance of key in .ini file and return them. If the
 * key is not specified, the proposal will be returned. If the proposal is -1,
 * but the key is not in the list, an error will be thrown.
 */
int countKeys(const char *key, const int proposal)
{
	/* make lowercase copy of the key */
	char *keyLower = malloc(strlen(key) + 1);
	strcpy(keyLower, key);
	for (int i = 0; i < strlen(keyLower); ++i) {
		keyLower[i] = tolower(keyLower[i]);
	}

	/* search for occourances of key */
	int count = 0;
	for (cmd_t *currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		if (!strcmp(keyLower, currCmd->key)) {
			++count;
		}
	}

	if (count == 0) {
		if (proposal == -1) {
			printf("Inifile missing necessary keyword: %s\n",
					keyLower);
			exit(1);
		} else {
			return proposal;
		}
	} else {
		return count;
	}
}

/*
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 */
int getInt(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	int value = strtol(valueStr, NULL, 10);
	printf("| %19s = %27d (%s)\n", key, value, defMsg);
	return value;
}

/*
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown.
 */
double getDbl(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	double value = strtod(valueStr, NULL);
	printf("| %19s = %27g (%s)\n", key, value, defMsg);
	return value;
}

/*
 * Find the key in the command list and return the corresponding value. If the
 * key is not specified, the proposal will be returned. If the proposal is
 * NULL, but the key is not in the list, an error will be thrown. The value in
 * the ini file is accepted as True, if it is a 'T', otherwise it is false.
 */
bool getBool(const char *key, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	bool value;
	if (!strcmp(valueStr, "T")) {
		value = true;
	} else {
		value = false;
	}
	printf("| %19s = %27s (%s)\n", key, (value ? "T" : "F"), defMsg);
	return value;
}

/*
 * Find the key in the command list and return the corresponding integer array.
 * If the key is not specified, the proposal will be returned. If the proposal
 * is NULL, but the key is not in the list, an error will be thrown. The value
 * in the ini file is accepted as True, if it is a 'T', otherwise it is false.
 */
int *getIntArray(const char *key, const int N, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	int *value = malloc(N * sizeof(int));

	char *tok = strtok(valueStr, ",");
	for (int i = 0; i < N; ++i) {
		value[i] = strtol(tok, NULL, 10);
		tok = strtok(NULL, ",");
	}

	//printf("| %17s[%d] = {%d", key, N, value[0]);
	//for (int i = 1; i < N; ++i) {
	//	printf(", %d", value[i]);
	//}
	//printf("} |%s|\n", defMsg);

	char line[STRLEN] = "";
	sprintf(line, "{%d", value[0]);
	for (int i = 1; i < N; ++i) {
		sprintf(line + strlen(line), ", %d", value[i]);
	}
	sprintf(line + strlen(line), "}");

	printf("| %18s[%d] = %27s (%s)\n", key, N, line, defMsg);
	return value;
}

/*
 * Find the key in the command list and return the corresponding double array.
 * If the key is not specified, the proposal will be returned. If the proposal
 * is NULL, but the key is not in the list, an error will be thrown. The value
 * in the ini file is accepted as True, if it is a 'T', otherwise it is false.
 */
double *getDblArray(const char *key, const int N, const char *proposal)
{
	char defMsg[8];
	char *valueStr = findCmd(key, defMsg, proposal);
	double *value = malloc(N * sizeof(double));

	char *tok = strtok(valueStr, ",");
	for (int i = 0; i < N; ++i) {
		value[i] = strtod(tok, NULL);
		tok = strtok(NULL, ",");
	}

	char line[STRLEN] = "";
	sprintf(line, "{%g", value[0]);
	for (int i = 1; i < N; ++i) {
		sprintf(line + strlen(line), ", %g", value[i]);
	}
	sprintf(line + strlen(line), "}");

	printf("| %18s[%d] = %27s (%s)\n", key, N, line, defMsg);
	return value;
}

/*
 * Savely delete all the elements in the list. The key and value strings need
 * be freeed seperately, because cmd_t only holds the pointers to those strings.
 */
void freeCmds(void)
{
	cmd_t *currCmd = firstCmd;
	do {
		currCmd = currCmd->next;
		free(currCmd->prev->key);
		free(currCmd->prev->value);
		free(currCmd->prev);
	} while (currCmd->next);
	free(currCmd->key);
	free(currCmd->value);
	free(currCmd);
}

/*
 * Print out all remaining commands in the list.
 */
void ignoredCmds(void)
{
	printf("The following commands were ignored:\n");
	for (cmd_t *currCmd = firstCmd; currCmd; currCmd = currCmd->next) {
		if (strlen(currCmd->value) > 0) {
			printf("\"%s\" = \"%s\"\n", currCmd->key, currCmd->value);
		} else {
			printf("[\"%s\" = \"%s\"]\n", currCmd->key, currCmd->value);
		}
	}

	/* Remove the remaining commands from memory. */
	freeCmds();
}
