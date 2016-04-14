#ifndef __UTEST_NB_CONTAINER_BOT_H__
#define __UTEST_NB_CONTAINER_BOT_H__

void cunit_nb_container_bot_array(void);
void cunit_nb_container_bot_QUEUE(void);
void cunit_nb_container_bot_STACK(void);
void cunit_nb_container_bot_SORTED(void);
void cunit_nb_container_bot_HASH(void);
void cunit_nb_container_bot_HEAP(void);
void cunit_nb_container_bot_iterator_QUEUE(void);
void cunit_nb_container_bot_iterator_STACK(void);
void cunit_nb_container_bot_iterator_SORTED(void);
void cunit_nb_container_bot_iterator_HASH(void);
void cunit_nb_container_bot_iterator_HEAP(void);

static void cunit_suites_nb_container_bot(void)
{
	cunit_nb_container_bot_array();
	cunit_nb_container_bot_QUEUE();
	cunit_nb_container_bot_STACK();
	cunit_nb_container_bot_SORTED();
	cunit_nb_container_bot_HASH();
	cunit_nb_container_bot_HEAP();
	cunit_nb_container_bot_iterator_QUEUE();
	cunit_nb_container_bot_iterator_STACK();
	cunit_nb_container_bot_iterator_SORTED();
	cunit_nb_container_bot_iterator_HASH();
	cunit_nb_container_bot_iterator_HEAP();
}

#endif
